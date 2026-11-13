/**
 * ALTD_WPhysics — 2D tunnel section of the AeroLab 3D fractional-step MAC scheme.
 *
 * Same ingredients as NavierStokes3D (tunnel boundaries, staggered faces,
 * density-weighted Poisson with p = 0 on the outlet face, smoothstep inlet
 * ramp from rest, air at 1.204 kg/m³ and 1.81e-5 Pa·s):
 *
 *   advect → diffuse → ∇²p = −ρ ∇·u* / Δt → u ← u* − (Δt/ρ) ∇p
 *
 * The browser integrator advects with semi-Lagrangian backtrace instead of
 * AeroLab's conservative upwind/MUSCL flux. First-order upwind on a phone-sized
 * mesh has numerical viscosity large enough to erase a Re 150 wake. Pressure is
 * relaxed with SOR, not the native PCG + geometric multigrid at 1e-7.
 * Coefficients are 2D, per unit depth, and are not the unverified 3D Cd.
 */

export type BodyId = "cylinder" | "cube" | "airfoil" | "funnel" | "torus" | "helix";

export type FieldMode = "smoke" | "speed" | "pressure" | "vorticity";

export interface TunnelSpec {
  body: BodyId;
  inlet: number;
  angleDeg: number;
  diameter: number;
  nx: number;
  ny: number;
}

export interface TunnelFrame {
  time: number;
  dt: number;
  steps: number;
  reynolds: number;
  inletNow: number;
  cd: number;
  cl: number;
  divergence: number;
  residual: number;
  maxSpeed: number;
  pressureMin: number;
  pressureMax: number;
  vorticityMax: number;
  warning: string;
}

export interface ScalarField {
  values: Float64Array;
  min: number;
  max: number;
  /** When set, t = 0.5 at this value (AeroLab pressure map). */
  reference: number | null;
  unit: string;
}

interface Particle {
  x: number;
  y: number;
  u: number;
}

const RHO = 1.204;
const MU = 1.81e-5;
const NU = MU / RHO;

export const FLUID = { density: RHO, viscosity: MU, nu: NU };

export function reynolds(inlet: number, diameter: number): number {
  return (Math.abs(inlet) * diameter) / NU;
}

function smoothstep(a: number): number {
  const x = a < 0 ? 0 : a > 1 ? 1 : a;
  return x * x * (3 - 2 * x);
}

export function containsBody(
  body: BodyId,
  x: number,
  y: number,
  cx: number,
  cy: number,
  angleRad: number,
  diameter: number,
): boolean {
  const c = Math.cos(angleRad);
  const s = Math.sin(angleRad);
  const dx = x - cx;
  const dy = y - cy;
  const xl = c * dx + s * dy;
  const yl = -s * dx + c * dy;
  const D = diameter;
  if (body === "cylinder") {
    const r = D * 0.5;
    return xl * xl + yl * yl <= r * r;
  }
  if (body === "cube") {
    const h = D * 0.5;
    return Math.abs(xl) <= h && Math.abs(yl) <= h;
  }
  if (body === "airfoil") {
    const chord = D;
    const xi = (xl + chord * 0.5) / chord;
    if (xi <= 0 || xi >= 1) return false;
    const t = 0.12;
    const rt = Math.sqrt(xi);
    const yt =
      5 *
      t *
      chord *
      (0.2969 * rt - 0.126 * xi - 0.3516 * xi * xi + 0.2843 * xi * xi * xi - 0.1015 * xi * xi * xi * xi);
    return Math.abs(yl) <= yt;
  }
  if (body === "funnel") {
    const L = D;
    const t = (xl + L * 0.5) / L;
    if (t < 0 || t > 1) return false;
    const half = L * (0.46 * (1 - t) + 0.14 * t);
    return Math.abs(yl) <= half;
  }
  if (body === "torus") {
    const r = D * 0.2;
    const sep = D * 0.3;
    const d1 = xl * xl + (yl - sep) * (yl - sep);
    const d2 = xl * xl + (yl + sep) * (yl + sep);
    return d1 <= r * r || d2 <= r * r;
  }
  const radius = D * 0.15;
  const r2 = radius * radius;
  for (let k = 0; k < 4; k += 1) {
    const ox = -1.5 * D * 0.4 + k * D * 0.4;
    const oy = (k % 2 === 0 ? 1 : -1) * D * 0.24;
    const ddx = xl - ox;
    const ddy = yl - oy;
    if (ddx * ddx + ddy * ddy <= r2) return true;
  }
  return false;
}

export class WindTunnel {
  nx = 120;
  ny = 60;
  lx = 0.2;
  ly = 0.1;
  dx = 0.2 / 120;
  dy = 0.1 / 60;

  body: BodyId = "cylinder";
  inlet = 0.113;
  angleRad = 0;
  diameter = 0.02;
  cx = 0.072;
  cy = 0.05;

  time = 0;
  steps = 0;
  ramp = 0.18;

  readonly solid: Uint8Array;
  readonly U: Float64Array;
  readonly V: Float64Array;
  readonly p: Float64Array;

  private Us: Float64Array;
  private Vs: Float64Array;
  private rhs: Float64Array;
  private scratch: Float64Array;
  private r: Float64Array;
  private z: Float64Array;
  private dir: Float64Array;
  private Ad: Float64Array;
  private diag: Float64Array;
  particles: Particle[] = [];
  warning = "";

  private emaCd = 0;
  private emaCl = 0;
  private emaReady = false;

  constructor(spec: TunnelSpec) {
    this.nx = spec.nx;
    this.ny = spec.ny;
    this.dx = this.lx / this.nx;
    this.dy = this.ly / this.ny;
    const nc = this.nx * this.ny;
    this.solid = new Uint8Array(nc);
    this.U = new Float64Array((this.nx + 1) * this.ny);
    this.V = new Float64Array(this.nx * (this.ny + 1));
    this.p = new Float64Array(nc);
    this.Us = new Float64Array(this.U.length);
    this.Vs = new Float64Array(this.V.length);
    this.rhs = new Float64Array(nc);
    this.scratch = new Float64Array(nc);
    this.r = new Float64Array(nc);
    this.z = new Float64Array(nc);
    this.dir = new Float64Array(nc);
    this.Ad = new Float64Array(nc);
    this.diag = new Float64Array(nc);
    this.rebuild(spec);
  }

  iu(i: number, j: number): number {
    return j * (this.nx + 1) + i;
  }

  iv(i: number, j: number): number {
    return j * this.nx + i;
  }

  ic(i: number, j: number): number {
    return j * this.nx + i;
  }

  cellFluid(i: number, j: number): boolean {
    if (i < 0 || j < 0 || i >= this.nx || j >= this.ny) return false;
    return this.solid[this.ic(i, j)] === 0;
  }

  openU(i: number, j: number): boolean {
    if (j < 0 || j >= this.ny || i < 0 || i > this.nx) return false;
    if (i === 0) return this.cellFluid(0, j);
    if (i === this.nx) return this.cellFluid(this.nx - 1, j);
    return this.cellFluid(i - 1, j) && this.cellFluid(i, j);
  }

  openV(i: number, j: number): boolean {
    if (i < 0 || i >= this.nx || j < 0 || j > this.ny) return false;
    if (j === 0 || j === this.ny) return false;
    return this.cellFluid(i, j - 1) && this.cellFluid(i, j);
  }

  /** Returns a placement warning, or "" when the mask is usable. */
  rebuild(spec: TunnelSpec): string {
    this.body = spec.body;
    this.inlet = spec.inlet;
    this.angleRad = (-spec.angleDeg * Math.PI) / 180;
    this.diameter = spec.diameter;
    this.cx = this.lx * 0.36;
    this.cy = this.ly * 0.5 + 0.08 * this.dy;
    this.ramp = this.diameter / Math.max(this.inlet, 1e-6);
    this.time = 0;
    this.steps = 0;
    this.emaReady = false;
    this.emaCd = 0;
    this.emaCl = 0;
    this.warning = "";
    this.U.fill(0);
    this.V.fill(0);
    this.p.fill(0);

    for (let j = 0; j < this.ny; j += 1) {
      for (let i = 0; i < this.nx; i += 1) {
        const x = (i + 0.5) * this.dx;
        const y = (j + 0.5) * this.dy;
        this.solid[this.ic(i, j)] = containsBody(
          this.body,
          x,
          y,
          this.cx,
          this.cy,
          this.angleRad,
          this.diameter,
        )
          ? 1
          : 0;
      }
    }

    const reach = new Uint8Array(this.nx * this.ny);
    const stack: number[] = [];
    for (let j = 0; j < this.ny; j += 1) {
      const id = this.ic(this.nx - 1, j);
      if (this.solid[id] === 0) {
        reach[id] = 1;
        stack.push(id);
      }
    }
    while (stack.length) {
      const a = stack.pop() as number;
      const i = a % this.nx;
      const j = (a / this.nx) | 0;
      const nbrs = [
        [i - 1, j],
        [i + 1, j],
        [i, j - 1],
        [i, j + 1],
      ];
      for (const [x, y] of nbrs) {
        if (!this.cellFluid(x, y)) continue;
        const b = this.ic(x, y);
        if (reach[b]) continue;
        reach[b] = 1;
        stack.push(b);
      }
    }
    let fluid = 0;
    for (let a = 0; a < this.solid.length; a += 1) {
      if (this.solid[a] === 0 && reach[a] === 0) this.solid[a] = 1;
      if (this.solid[a] === 0) fluid += 1;
    }

    let blocked = false;
    for (let j = 0; j < this.ny; j += 1) {
      if (!this.cellFluid(0, j)) blocked = true;
    }
    if (blocked) this.warning = "The body touches the inlet. Reduce size or angle.";
    else if (fluid < (this.nx * this.ny) / 8) this.warning = "The body blocks most of the section.";

    this.boundary(this.U, this.V);
    this.buildDiagonal();
    this.seedParticles();
    return this.warning;
  }

  setInlet(speed: number): void {
    this.inlet = speed;
    this.ramp = this.diameter / Math.max(speed, 1e-6);
  }

  imposed(atTime: number): number {
    const a = smoothstep(atTime / Math.max(this.ramp, 1e-9));
    return this.inlet * a;
  }

  private boundary(u: Float64Array, v: Float64Array): void {
    const uin = this.imposed(this.time);
    for (let j = 0; j < this.ny; j += 1) {
      for (let i = 0; i <= this.nx; i += 1) {
        const id = this.iu(i, j);
        if (!this.openU(i, j)) u[id] = 0;
        else if (i === 0) u[id] = uin;
      }
    }
    for (let j = 0; j <= this.ny; j += 1) {
      for (let i = 0; i < this.nx; i += 1) {
        if (!this.openV(i, j)) v[this.iv(i, j)] = 0;
      }
    }
  }

  private sampleU(x: number, y: number): number {
    if (x < 0) return this.imposed(this.time);
    const maxX = this.lx;
    const xx = x > maxX ? maxX : x;
    let yy = y;
    const y0 = 0.5 * this.dy;
    const y1 = this.ly - 0.5 * this.dy;
    if (yy < y0) yy = y0;
    if (yy > y1) yy = y1;
    let gx = xx / this.dx;
    let gy = yy / this.dy - 0.5;
    if (gx < 0) gx = 0;
    if (gx > this.nx) gx = this.nx;
    if (gy < 0) gy = 0;
    if (gy > this.ny - 1) gy = this.ny - 1;
    const i0 = Math.floor(gx);
    const j0 = Math.floor(gy);
    const i1 = i0 + 1 > this.nx ? this.nx : i0 + 1;
    const j1 = j0 + 1 > this.ny - 1 ? this.ny - 1 : j0 + 1;
    const tx = gx - i0;
    const ty = gy - j0;
    const a = this.U[this.iu(i0, j0)];
    const b = this.U[this.iu(i1, j0)];
    const c = this.U[this.iu(i0, j1)];
    const d = this.U[this.iu(i1, j1)];
    return (a * (1 - tx) + b * tx) * (1 - ty) + (c * (1 - tx) + d * tx) * ty;
  }

  private sampleV(x: number, y: number): number {
    if (y <= 0 || y >= this.ly) return 0;
    let xx = x;
    if (xx < 0.5 * this.dx) xx = 0.5 * this.dx;
    if (xx > this.lx - 0.5 * this.dx) xx = this.lx - 0.5 * this.dx;
    let gx = xx / this.dx - 0.5;
    let gy = y / this.dy;
    if (gx < 0) gx = 0;
    if (gx > this.nx - 1) gx = this.nx - 1;
    if (gy < 0) gy = 0;
    if (gy > this.ny) gy = this.ny;
    const i0 = Math.floor(gx);
    const j0 = Math.floor(gy);
    const i1 = i0 + 1 > this.nx - 1 ? this.nx - 1 : i0 + 1;
    const j1 = j0 + 1 > this.ny ? this.ny : j0 + 1;
    const tx = gx - i0;
    const ty = gy - j0;
    const a = this.V[this.iv(i0, j0)];
    const b = this.V[this.iv(i1, j0)];
    const c = this.V[this.iv(i0, j1)];
    const d = this.V[this.iv(i1, j1)];
    return (a * (1 - tx) + b * tx) * (1 - ty) + (c * (1 - tx) + d * tx) * ty;
  }

  velocityAt(x: number, y: number): { u: number; v: number } {
    if (x < 0 || y < 0 || x > this.lx || y > this.ly) return { u: 0, v: 0 };
    if (containsBody(this.body, x, y, this.cx, this.cy, this.angleRad, this.diameter)) {
      return { u: 0, v: 0 };
    }
    return { u: this.sampleU(x, y), v: this.sampleV(x, y) };
  }

  private chooseDt(): number {
    let umax = Math.abs(this.inlet);
    let vmax = 0;
    for (let i = 0; i < this.U.length; i += 1) {
      const a = Math.abs(this.U[i]);
      if (a > umax) umax = a;
    }
    for (let i = 0; i < this.V.length; i += 1) {
      const a = Math.abs(this.V[i]);
      if (a > vmax) vmax = a;
    }
    const rate =
      umax / this.dx +
      vmax / this.dy +
      2 * NU * (1 / (this.dx * this.dx) + 1 / (this.dy * this.dy));
    const dt = Math.min(0.02, 0.85 / Math.max(rate, 1e-9));
    return dt > 0 && Number.isFinite(dt) ? dt : 1e-4;
  }

  private advectDiffuse(dt: number): void {
    const uin = this.imposed(this.time);
    for (let j = 0; j < this.ny; j += 1) {
      for (let i = 0; i <= this.nx; i += 1) {
        const id = this.iu(i, j);
        if (!this.openU(i, j)) {
          this.Us[id] = 0;
          continue;
        }
        if (i === 0) {
          this.Us[id] = uin;
          continue;
        }
        const x = i * this.dx;
        const y = (j + 0.5) * this.dy;
        const uu = this.U[id];
        const vv = this.sampleV(x, y);
        this.Us[id] = this.sampleU(x - dt * uu, y - dt * vv);
      }
    }
    for (let j = 0; j <= this.ny; j += 1) {
      for (let i = 0; i < this.nx; i += 1) {
        const id = this.iv(i, j);
        if (!this.openV(i, j)) {
          this.Vs[id] = 0;
          continue;
        }
        const x = (i + 0.5) * this.dx;
        const y = j * this.dy;
        const uu = this.sampleU(x, y);
        const vv = this.V[id];
        this.Vs[id] = this.sampleV(x - dt * uu, y - dt * vv);
      }
    }

    const dx2 = this.dx * this.dx;
    const dy2 = this.dy * this.dy;
    const vis = NU * dt;
    for (let j = 0; j < this.ny; j += 1) {
      for (let i = 1; i < this.nx; i += 1) {
        if (!this.openU(i, j)) continue;
        const id = this.iu(i, j);
        const c = this.Us[id];
        const e = this.openU(i + 1, j) ? this.Us[this.iu(i + 1, j)] : c;
        const w = this.openU(i - 1, j) ? this.Us[this.iu(i - 1, j)] : 0;
        const n = j + 1 < this.ny && this.openU(i, j + 1) ? this.Us[this.iu(i, j + 1)] : c;
        const s = j - 1 >= 0 && this.openU(i, j - 1) ? this.Us[this.iu(i, j - 1)] : c;
        const lap = (e - 2 * c + w) / dx2 + (n - 2 * c + s) / dy2;
        this.U[id] = c + vis * lap;
      }
    }
    for (let j = 1; j < this.ny; j += 1) {
      for (let i = 0; i < this.nx; i += 1) {
        if (!this.openV(i, j)) continue;
        const id = this.iv(i, j);
        const c = this.Vs[id];
        const e = i + 1 < this.nx && this.openV(i + 1, j) ? this.Vs[this.iv(i + 1, j)] : c;
        const w = i - 1 >= 0 && this.openV(i - 1, j) ? this.Vs[this.iv(i - 1, j)] : c;
        const n = this.openV(i, j + 1) ? this.Vs[this.iv(i, j + 1)] : 0;
        const s = this.openV(i, j - 1) ? this.Vs[this.iv(i, j - 1)] : 0;
        const lap = (e - 2 * c + w) / dx2 + (n - 2 * c + s) / dy2;
        this.V[id] = c + vis * lap;
      }
    }
    this.Us.set(this.U);
    this.Vs.set(this.V);
    for (let j = 0; j < this.ny; j += 1) {
      if (!this.openU(this.nx, j)) continue;
      this.Us[this.iu(this.nx, j)] = this.Us[this.iu(this.nx - 1, j)];
    }
    this.boundary(this.Us, this.Vs);
  }

  private formRhs(dt: number): void {
    const ix = 1 / this.dx;
    const iy = 1 / this.dy;
    for (let j = 0; j < this.ny; j += 1) {
      for (let i = 0; i < this.nx; i += 1) {
        const id = this.ic(i, j);
        if (this.solid[id]) {
          this.rhs[id] = 0;
          continue;
        }
        const div =
          (this.Us[this.iu(i + 1, j)] - this.Us[this.iu(i, j)]) * ix +
          (this.Vs[this.iv(i, j + 1)] - this.Vs[this.iv(i, j)]) * iy;
        this.rhs[id] = (-RHO * div) / dt;
      }
    }
  }

  private buildDiagonal(): void {
    const ax = 1 / (this.dx * this.dx);
    const ay = 1 / (this.dy * this.dy);
    for (let j = 0; j < this.ny; j += 1) {
      for (let i = 0; i < this.nx; i += 1) {
        const id = this.ic(i, j);
        if (this.solid[id]) {
          this.diag[id] = 1;
          this.p[id] = 0;
          continue;
        }
        let diag = 0;
        if (i > 0 && this.cellFluid(i - 1, j)) diag += ax;
        if (i < this.nx - 1 && this.cellFluid(i + 1, j)) diag += ax;
        if (i === this.nx - 1) diag += 2 * ax;
        if (j > 0 && this.cellFluid(i, j - 1)) diag += ay;
        if (j < this.ny - 1 && this.cellFluid(i, j + 1)) diag += ay;
        this.diag[id] = diag > 0 ? diag : 1;
      }
    }
  }

  private applyA(x: Float64Array, out: Float64Array): void {
    const ax = 1 / (this.dx * this.dx);
    const ay = 1 / (this.dy * this.dy);
    const nx = this.nx;
    const ny = this.ny;
    for (let j = 0; j < ny; j += 1) {
      for (let i = 0; i < nx; i += 1) {
        const id = this.ic(i, j);
        if (this.solid[id]) {
          out[id] = 0;
          continue;
        }
        const val = x[id];
        let result = 0;
        if (i > 0 && this.cellFluid(i - 1, j)) result += ax * (val - x[this.ic(i - 1, j)]);
        if (i < nx - 1 && this.cellFluid(i + 1, j)) result += ax * (val - x[this.ic(i + 1, j)]);
        if (i === nx - 1) result += 2 * ax * val;
        if (j > 0 && this.cellFluid(i, j - 1)) result += ay * (val - x[this.ic(i, j - 1)]);
        if (j < ny - 1 && this.cellFluid(i, j + 1)) result += ay * (val - x[this.ic(i, j + 1)]);
        out[id] = result;
      }
    }
  }

  /** Jacobi-preconditioned CG, capped. Same operator as AeroLab's tunnel Laplacian. */
  private pressureSolve(maxIter: number): number {
    this.applyA(this.p, this.Ad);
    let r2 = 0;
    let b2 = 0;
    let rz = 0;
    for (let a = 0; a < this.p.length; a += 1) {
      if (this.solid[a]) {
        this.r[a] = 0;
        this.z[a] = 0;
        this.dir[a] = 0;
        this.p[a] = 0;
        continue;
      }
      const res = this.rhs[a] - this.Ad[a];
      this.r[a] = res;
      this.z[a] = res / this.diag[a];
      this.dir[a] = this.z[a];
      r2 += res * res;
      b2 += this.rhs[a] * this.rhs[a];
      rz += res * this.z[a];
    }
    const bnorm = Math.sqrt(b2);
    const threshold = Math.max(1e-8, bnorm * 1e-4);
    if (Math.sqrt(r2) <= threshold) return Math.sqrt(r2) / Math.max(1e-14, bnorm);

    let residual = Math.sqrt(r2) / Math.max(1e-14, bnorm);
    for (let it = 0; it < maxIter; it += 1) {
      this.applyA(this.dir, this.Ad);
      let pAp = 0;
      for (let a = 0; a < this.p.length; a += 1) pAp += this.dir[a] * this.Ad[a];
      if (!(pAp > 0) || !Number.isFinite(pAp) || !Number.isFinite(rz)) break;
      const alpha = rz / pAp;
      r2 = 0;
      for (let a = 0; a < this.p.length; a += 1) {
        if (this.solid[a]) continue;
        this.p[a] += alpha * this.dir[a];
        this.r[a] -= alpha * this.Ad[a];
        r2 += this.r[a] * this.r[a];
      }
      residual = Math.sqrt(r2) / Math.max(1e-14, bnorm);
      if (Math.sqrt(r2) <= threshold) break;
      let nextRz = 0;
      for (let a = 0; a < this.p.length; a += 1) {
        if (this.solid[a]) {
          this.z[a] = 0;
          continue;
        }
        this.z[a] = this.r[a] / this.diag[a];
        nextRz += this.r[a] * this.z[a];
      }
      const beta = nextRz / rz;
      for (let a = 0; a < this.p.length; a += 1) {
        if (this.solid[a]) {
          this.dir[a] = 0;
          continue;
        }
        this.dir[a] = this.z[a] + beta * this.dir[a];
      }
      rz = nextRz;
    }
    return residual;
  }

  private project(dt: number): void {
    const scale = dt / RHO;
    for (let j = 0; j < this.ny; j += 1) {
      for (let i = 1; i <= this.nx; i += 1) {
        if (!this.openU(i, j)) continue;
        const grad =
          i === this.nx
            ? (-this.p[this.ic(i - 1, j)] * 2) / this.dx
            : (this.p[this.ic(i, j)] - this.p[this.ic(i - 1, j)]) / this.dx;
        this.Us[this.iu(i, j)] -= scale * grad;
      }
    }
    for (let j = 1; j < this.ny; j += 1) {
      for (let i = 0; i < this.nx; i += 1) {
        if (!this.openV(i, j)) continue;
        const grad = (this.p[this.ic(i, j)] - this.p[this.ic(i, j - 1)]) / this.dy;
        this.Vs[this.iv(i, j)] -= scale * grad;
      }
    }
    this.boundary(this.Us, this.Vs);
    this.U.set(this.Us);
    this.V.set(this.Vs);
  }

  private forces(): { cd: number; cl: number } {
    let fx = 0;
    let fy = 0;
    const dAx = this.dy;
    const dAy = this.dx;
    for (let j = 0; j < this.ny; j += 1) {
      for (let i = 1; i < this.nx; i += 1) {
        const left = this.solid[this.ic(i - 1, j)];
        const right = this.solid[this.ic(i, j)];
        if (left === right) continue;
        const fluidIsRight = right === 0;
        const fluidI = fluidIsRight ? i : i - 1;
        const n = fluidIsRight ? 1 : -1;
        const pf = this.p[this.ic(fluidI, j)];
        fx += -pf * n * dAx;
        const vt = 0.5 * (this.V[this.iv(fluidI, j)] + this.V[this.iv(fluidI, j + 1)]);
        fy += ((MU * vt) / (0.5 * this.dx)) * dAx;
      }
    }
    for (let j = 1; j < this.ny; j += 1) {
      for (let i = 0; i < this.nx; i += 1) {
        const below = this.solid[this.ic(i, j - 1)];
        const above = this.solid[this.ic(i, j)];
        if (below === above) continue;
        const fluidIsAbove = above === 0;
        const fluidJ = fluidIsAbove ? j : j - 1;
        const n = fluidIsAbove ? 1 : -1;
        const pf = this.p[this.ic(i, fluidJ)];
        fy += -pf * n * dAy;
        const ut = 0.5 * (this.U[this.iu(i, fluidJ)] + this.U[this.iu(i + 1, fluidJ)]);
        const shear = MU * (ut / (0.5 * this.dy));
        fx += shear * dAy;
      }
    }
    const q = 0.5 * RHO * this.inlet * this.inlet * Math.max(this.diameter, 1e-6);
    const cd = q > 0 ? fx / q : 0;
    const cl = q > 0 ? fy / q : 0;
    if (!this.emaReady) {
      this.emaCd = cd;
      this.emaCl = cl;
      this.emaReady = true;
    } else {
      this.emaCd = this.emaCd * 0.92 + cd * 0.08;
      this.emaCl = this.emaCl * 0.92 + cl * 0.08;
    }
    return { cd: this.emaCd, cl: this.emaCl };
  }

  private divergence(): number {
    let sum = 0;
    let n = 0;
    const ix = 1 / this.dx;
    const iy = 1 / this.dy;
    for (let j = 0; j < this.ny; j += 1) {
      for (let i = 0; i < this.nx; i += 1) {
        if (!this.cellFluid(i, j)) continue;
        const div =
          (this.U[this.iu(i + 1, j)] - this.U[this.iu(i, j)]) * ix +
          (this.V[this.iv(i, j + 1)] - this.V[this.iv(i, j)]) * iy;
        sum += div * div;
        n += 1;
      }
    }
    return Math.sqrt(sum / Math.max(1, n));
  }

  private seedParticles(): void {
    const count = Math.min(1100, Math.max(420, Math.floor(this.nx * this.ny * 0.12)));
    this.particles = new Array(count);
    for (let i = 0; i < count; i += 1) this.particles[i] = this.spawn(i);
  }

  private spawn(salt: number): Particle {
    const y = this.dy * (1 + ((salt * 97) % Math.max(1, this.ny - 2)));
    const x = this.dx * (1.2 + (salt % 5) * 0.7);
    return { x, y, u: this.imposed(this.time) };
  }

  private advectParticles(dt: number): void {
    const parts = this.particles;
    for (let i = 0; i < parts.length; i += 1) {
      const p = parts[i];
      const a = this.velocityAt(p.x, p.y);
      const mx = p.x + 0.5 * dt * a.u;
      const my = p.y + 0.5 * dt * a.v;
      const b = this.velocityAt(mx, my);
      p.x += dt * b.u;
      p.y += dt * b.v;
      p.u = b.u;
      const inside = containsBody(this.body, p.x, p.y, this.cx, this.cy, this.angleRad, this.diameter);
      if (inside || p.x < 0 || p.x > this.lx || p.y < this.dy || p.y > this.ly - this.dy) {
        parts[i] = this.spawn(i + this.steps);
      }
    }
  }

  step(sorIters = 28): TunnelFrame {
    if (this.warning) {
      return this.frame(0, 0, 0, 0);
    }
    const dt = this.chooseDt();
    this.advectDiffuse(dt);
    this.formRhs(dt);
    const residual = this.pressureSolve(sorIters);
    this.project(dt);
    this.time += dt;
    this.steps += 1;
    this.advectParticles(dt);
    const div = this.divergence();
    const { cd, cl } = this.forces();
    return this.frame(dt, residual, div, cd, cl);
  }

  private frame(dt: number, residual: number, divergence: number, cd = 0, cl = 0): TunnelFrame {
    let maxSpeed = 0;
    let pMin = Infinity;
    let pMax = -Infinity;
    let wMax = 0;
    for (let j = 0; j < this.ny; j += 1) {
      for (let i = 0; i < this.nx; i += 1) {
        if (!this.cellFluid(i, j)) continue;
        const u = 0.5 * (this.U[this.iu(i, j)] + this.U[this.iu(i + 1, j)]);
        const v = 0.5 * (this.V[this.iv(i, j)] + this.V[this.iv(i, j + 1)]);
        const sp = Math.hypot(u, v);
        if (sp > maxSpeed) maxSpeed = sp;
        const pv = this.p[this.ic(i, j)];
        if (pv < pMin) pMin = pv;
        if (pv > pMax) pMax = pv;
        const vort =
          (this.V[this.iv(i + 1 > this.nx - 1 ? this.nx - 1 : i + 1, j)] -
            this.V[this.iv(i, j)]) /
            this.dx -
          (this.U[this.iu(i, j + 1 > this.ny - 1 ? this.ny - 1 : j + 1)] - this.U[this.iu(i, j)]) /
            this.dy;
        const aw = Math.abs(vort);
        if (aw > wMax) wMax = aw;
      }
    }
    if (!Number.isFinite(pMin)) {
      pMin = 0;
      pMax = 0;
    }
    return {
      time: this.time,
      dt,
      steps: this.steps,
      reynolds: reynolds(this.inlet, this.diameter),
      inletNow: this.imposed(this.time),
      cd,
      cl,
      divergence,
      residual,
      maxSpeed,
      pressureMin: pMin,
      pressureMax: pMax,
      vorticityMax: wMax,
      warning: this.warning,
    };
  }

  fillScalar(mode: Exclude<FieldMode, "smoke">): ScalarField {
    const values = this.scratch;
    let min = Infinity;
    let max = -Infinity;
    for (let j = 0; j < this.ny; j += 1) {
      for (let i = 0; i < this.nx; i += 1) {
        const id = this.ic(i, j);
        if (!this.cellFluid(i, j)) {
          values[id] = 0;
          continue;
        }
        const u = 0.5 * (this.U[this.iu(i, j)] + this.U[this.iu(i + 1, j)]);
        const v = 0.5 * (this.V[this.iv(i, j)] + this.V[this.iv(i, j + 1)]);
        let s = 0;
        if (mode === "speed") s = Math.hypot(u, v);
        else if (mode === "pressure") s = this.p[id];
        else {
          const ip = i + 1 > this.nx - 1 ? this.nx - 1 : i + 1;
          const jp = j + 1 > this.ny - 1 ? this.ny - 1 : j + 1;
          s = (this.V[this.iv(ip, j)] - this.V[this.iv(i, j)]) / this.dx - (this.U[this.iu(i, jp)] - this.U[this.iu(i, j)]) / this.dy;
        }
        values[id] = s;
        if (s < min) min = s;
        if (s > max) max = s;
      }
    }
    if (!Number.isFinite(min) || min === max) {
      min = 0;
      max = 1;
    }
    if (mode === "pressure" || mode === "vorticity") {
      const ext = Math.max(Math.abs(min), Math.abs(max), 1e-12);
      return { values, min: -ext, max: ext, reference: 0, unit: mode === "pressure" ? "Pa" : "1/s" };
    }
    return { values, min, max, reference: null, unit: "m/s" };
  }
}
