import { useEffect, useRef, useState } from "react";
import { Pause, Play, RotateCcw } from "lucide-react";
import { buildLut, type MapKind } from "@/lib/wind/colormap";
import {
  FLUID,
  WindTunnel,
  containsBody,
  type BodyId,
  type FieldMode,
  type TunnelFrame,
  type TunnelSpec,
} from "@/lib/wind/solver";

const BODIES: { id: BodyId; label: string }[] = [
  { id: "cylinder", label: "Cylinder" },
  { id: "cube", label: "Cube" },
  { id: "airfoil", label: "Airfoil" },
  { id: "funnel", label: "Funnel" },
  { id: "torus", label: "Torus" },
  { id: "helix", label: "Helix" },
];

const MODES: { id: FieldMode; label: string }[] = [
  { id: "smoke", label: "Smoke" },
  { id: "speed", label: "Speed" },
  { id: "pressure", label: "Pressure" },
  { id: "vorticity", label: "Vorticity" },
];

interface Preset {
  id: string;
  label: string;
  body: BodyId;
  inlet: number;
  angle: number;
  diameter: number;
}

const PRESETS: Preset[] = [
  { id: "re150", label: "Re 150 cylinder", body: "cylinder", inlet: 0.113, angle: 0, diameter: 0.02 },
  { id: "cube", label: "Cube", body: "cube", inlet: 0.15, angle: 8, diameter: 0.02 },
  { id: "foil", label: "Airfoil 8°", body: "airfoil", inlet: 0.22, angle: 8, diameter: 0.03 },
  { id: "funnel", label: "Funnel", body: "funnel", inlet: 0.16, angle: 0, diameter: 0.032 },
  { id: "gap", label: "Torus gap", body: "torus", inlet: 0.18, angle: 0, diameter: 0.036 },
  { id: "helix", label: "Helix", body: "helix", inlet: 0.16, angle: 10, diameter: 0.03 },
];

const EMPTY: TunnelFrame = {
  time: 0,
  dt: 0,
  steps: 0,
  reynolds: 150,
  inletNow: 0,
  cd: 0,
  cl: 0,
  divergence: 0,
  residual: 1,
  maxSpeed: 0,
  pressureMin: 0,
  pressureMax: 0,
  vorticityMax: 0,
  warning: "",
};

function gridForWidth(): { nx: number; ny: number } {
  return { nx: 96, ny: 48 };
}

function fmt(n: number, digits: number): string {
  if (!Number.isFinite(n)) return "—";
  const a = Math.abs(n);
  if (a !== 0 && (a < 0.001 || a >= 10000)) return n.toExponential(2);
  return n.toFixed(digits);
}

function shade(value: number, min: number, max: number, reference: number | null): number {
  let t: number;
  if (reference !== null && reference >= min && reference <= max) {
    if (value <= reference) t = reference > min ? 0.5 * ((value - min) / (reference - min)) : 0.5;
    else t = max > reference ? 0.5 + 0.5 * ((value - reference) / (max - reference)) : 0.5;
  } else {
    t = max > min ? (value - min) / (max - min) : 0.5;
  }
  if (t < 0) t = 0;
  if (t > 1) t = 1;
  return (t * 255) | 0;
}

function worldPoint(
  xl: number,
  yl: number,
  cx: number,
  cy: number,
  ang: number,
): { x: number; y: number } {
  const c = Math.cos(ang);
  const s = Math.sin(ang);
  return { x: cx + c * xl - s * yl, y: cy + s * xl + c * yl };
}

function traceBody(ctx: CanvasRenderingContext2D, tunnel: WindTunnel, w: number, h: number): void {
  const { body, cx, cy, angleRad, diameter, lx, ly } = tunnel;
  const X = (x: number) => (x / lx) * w;
  const Y = (y: number) => h - (y / ly) * h;
  ctx.beginPath();
  if (body === "cylinder") {
    const steps = 48;
    for (let i = 0; i <= steps; i += 1) {
      const a = (i / steps) * Math.PI * 2;
      const p = worldPoint(Math.cos(a) * diameter * 0.5, Math.sin(a) * diameter * 0.5, cx, cy, angleRad);
      if (i === 0) ctx.moveTo(X(p.x), Y(p.y));
      else ctx.lineTo(X(p.x), Y(p.y));
    }
  } else if (body === "cube") {
    const h2 = diameter * 0.5;
    const corners = [
      [-h2, -h2],
      [h2, -h2],
      [h2, h2],
      [-h2, h2],
    ];
    corners.forEach(([xl, yl], i) => {
      const p = worldPoint(xl, yl, cx, cy, angleRad);
      if (i === 0) ctx.moveTo(X(p.x), Y(p.y));
      else ctx.lineTo(X(p.x), Y(p.y));
    });
    ctx.closePath();
  } else if (body === "airfoil") {
    const chord = diameter;
    const n = 36;
    const yt = (xi: number) => {
      const t = 0.12;
      const rt = Math.sqrt(Math.max(xi, 0));
      return (
        5 *
        t *
        chord *
        (0.2969 * rt - 0.126 * xi - 0.3516 * xi * xi + 0.2843 * xi * xi * xi - 0.1015 * xi * xi * xi * xi)
      );
    };
    for (let i = 0; i <= n; i += 1) {
      const xi = i / n;
      const p = worldPoint(-chord * 0.5 + xi * chord, yt(xi), cx, cy, angleRad);
      if (i === 0) ctx.moveTo(X(p.x), Y(p.y));
      else ctx.lineTo(X(p.x), Y(p.y));
    }
    for (let i = n; i >= 0; i -= 1) {
      const xi = i / n;
      const p = worldPoint(-chord * 0.5 + xi * chord, -yt(xi), cx, cy, angleRad);
      ctx.lineTo(X(p.x), Y(p.y));
    }
    ctx.closePath();
  } else if (body === "funnel") {
    const L = diameter;
    const n = 16;
    for (let i = 0; i <= n; i += 1) {
      const t = i / n;
      const xl = -L * 0.5 + t * L;
      const half = L * (0.46 * (1 - t) + 0.14 * t);
      const p = worldPoint(xl, half, cx, cy, angleRad);
      if (i === 0) ctx.moveTo(X(p.x), Y(p.y));
      else ctx.lineTo(X(p.x), Y(p.y));
    }
    for (let i = n; i >= 0; i -= 1) {
      const t = i / n;
      const xl = -L * 0.5 + t * L;
      const half = L * (0.46 * (1 - t) + 0.14 * t);
      const p = worldPoint(xl, -half, cx, cy, angleRad);
      ctx.lineTo(X(p.x), Y(p.y));
    }
    ctx.closePath();
  } else if (body === "torus") {
    const r = diameter * 0.2;
    const sep = diameter * 0.3;
    const steps = 28;
    for (const oy of [sep, -sep]) {
      for (let i = 0; i <= steps; i += 1) {
        const a = (i / steps) * Math.PI * 2;
        const p = worldPoint(Math.cos(a) * r, oy + Math.sin(a) * r, cx, cy, angleRad);
        if (i === 0) ctx.moveTo(X(p.x), Y(p.y));
        else ctx.lineTo(X(p.x), Y(p.y));
      }
      ctx.closePath();
    }
  } else {
    const radius = diameter * 0.15;
    const steps = 20;
    for (let k = 0; k < 4; k += 1) {
      const ox = -1.5 * diameter * 0.4 + k * diameter * 0.4;
      const oy = (k % 2 === 0 ? 1 : -1) * diameter * 0.24;
      for (let i = 0; i <= steps; i += 1) {
        const a = (i / steps) * Math.PI * 2;
        const p = worldPoint(ox + Math.cos(a) * radius, oy + Math.sin(a) * radius, cx, cy, angleRad);
        if (i === 0) ctx.moveTo(X(p.x), Y(p.y));
        else ctx.lineTo(X(p.x), Y(p.y));
      }
      ctx.closePath();
    }
  }
  ctx.fillStyle = "#0c161f";
  ctx.strokeStyle = "#d7e6ef";
  ctx.lineWidth = Math.max(1, w / 420);
  ctx.fill();
  ctx.stroke();
}

function drawStreamlines(ctx: CanvasRenderingContext2D, tunnel: WindTunnel, w: number, h: number): void {
  const seeds = 12;
  ctx.beginPath();
  for (let s = 0; s < seeds; s += 1) {
    let x = tunnel.dx * 2;
    let y = tunnel.dy + ((tunnel.ly - 2 * tunnel.dy) * s) / (seeds - 1);
    let moved = false;
    for (let k = 0; k < 48; k += 1) {
      const vel = tunnel.velocityAt(x, y);
      const sp = Math.hypot(vel.u, vel.v);
      if (sp < 1e-5) break;
      const ds = tunnel.dx * 0.85;
      x += (vel.u / sp) * ds;
      y += (vel.v / sp) * ds;
      if (x < 0 || y < 0 || x > tunnel.lx || y > tunnel.ly) break;
      if (containsBody(tunnel.body, x, y, tunnel.cx, tunnel.cy, tunnel.angleRad, tunnel.diameter)) break;
      const px = (x / tunnel.lx) * w;
      const py = h - (y / tunnel.ly) * h;
      if (!moved) {
        ctx.moveTo(px, py);
        moved = true;
      } else ctx.lineTo(px, py);
    }
  }
  ctx.strokeStyle = "rgba(215,230,239,0.55)";
  ctx.lineWidth = 1;
  ctx.stroke();
}

export function WindLab() {
  const canvasRef = useRef<HTMLCanvasElement | null>(null);
  const fieldRef = useRef<HTMLCanvasElement | null>(null);
  const tunnelRef = useRef<WindTunnel | null>(null);
  const pausedRef = useRef(false);
  const modeRef = useRef<FieldMode>("smoke");
  const luts = useRef({
    contrast: buildLut("contrast"),
    viridis: buildLut("viridis"),
    pressure: buildLut("pressure"),
  });

  const [spec, setSpec] = useState<TunnelSpec>({
    body: "cylinder",
    inlet: 0.113,
    angleDeg: 0,
    diameter: 0.02,
    nx: 96,
    ny: 48,
  });
  const [mode, setMode] = useState<FieldMode>("smoke");
  const [paused, setPaused] = useState(false);
  const [stats, setStats] = useState<TunnelFrame>(EMPTY);
  const [noteOpen, setNoteOpen] = useState(false);

  useEffect(() => {
    pausedRef.current = paused;
  }, [paused]);
  useEffect(() => {
    modeRef.current = mode;
  }, [mode]);

  useEffect(() => {
    const next = { ...spec, ...gridForWidth() };
    const tunnel = new WindTunnel(next);
    tunnelRef.current = tunnel;
    setSpec(next);
    setStats(tunnel.step(12));

    let raf = 0;
    let lastHud = 0;
    const loop = (now: number) => {
      const live = tunnelRef.current;
      const canvas = canvasRef.current;
      if (live && canvas) {
        if (!pausedRef.current && !live.warning) {
          const frame = live.step(40);
          if (now - lastHud > 120) {
            lastHud = now;
            setStats(frame);
          }
        }
        paint(canvas, live, modeRef.current, luts.current, fieldRef);
      }
      raf = requestAnimationFrame(loop);
    };
    raf = requestAnimationFrame(loop);
    return () => cancelAnimationFrame(raf);
    // Mount only — later spec edits go through the handlers so the loop is not reset.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  function applyGeometry(partial: Partial<TunnelSpec>) {
    const live = tunnelRef.current;
    if (!live) return;
    const next: TunnelSpec = {
      body: partial.body ?? spec.body,
      inlet: partial.inlet ?? spec.inlet,
      angleDeg: partial.angleDeg ?? spec.angleDeg,
      diameter: partial.diameter ?? spec.diameter,
      nx: live.nx,
      ny: live.ny,
    };
    const warning = live.rebuild(next);
    setSpec(next);
    setStats({ ...live.step(12), warning });
  }

  function applyInlet(inlet: number) {
    const live = tunnelRef.current;
    if (!live) return;
    live.setInlet(inlet);
    setSpec((prev) => ({ ...prev, inlet }));
  }

  function reset() {
    applyGeometry({});
  }

  const legendKind: MapKind = mode === "pressure" || mode === "vorticity" ? "pressure" : mode === "speed" ? "viridis" : "contrast";
  const legendLabel =
    mode === "pressure" ? "gauge pressure" : mode === "vorticity" ? "vorticity" : mode === "speed" ? "speed" : "streamwise smoke";

  return (
    <div className="mx-auto flex w-full max-w-5xl flex-col gap-4 px-4 py-4 sm:px-6 sm:py-6">
      <header className="flex flex-wrap items-end justify-between gap-3">
        <div>
          <p className="font-mono text-xs tracking-widest text-primary uppercase">Section solver</p>
          <h1 className="text-2xl font-semibold text-fg sm:text-3xl">ALTD_WPhysics</h1>
          <p className="max-w-xl text-sm text-pretty text-muted">
            A 2D wind-tunnel slice of the AeroLab 3D MAC scheme. Air flows left to right. The inlet ramps
            from rest.
          </p>
        </div>
        <button
          type="button"
          onClick={() => setPaused((v) => !v)}
          className="inline-flex min-h-11 items-center gap-2 rounded-lg border border-border bg-surface px-3 text-sm text-fg"
        >
          {paused ? <Play className="size-4" /> : <Pause className="size-4" />}
          {paused ? "Run" : "Pause"}
        </button>
      </header>

      <div className="overflow-hidden rounded-lg border border-border bg-surface">
        <div className="tunnel-frame relative">
          <canvas ref={canvasRef} className="absolute inset-0 h-full w-full" aria-label="Wind tunnel section" />
        </div>
        <dl className="grid grid-cols-2 gap-px border-t border-border bg-border sm:grid-cols-4">
          <Stat label="Re" value={fmt(stats.reynolds, 0)} hint="U D / ν" />
          <Stat label="U inlet" value={`${fmt(stats.inletNow, 3)} m/s`} hint={`target ${fmt(spec.inlet, 3)}`} />
          <Stat label="Cd" value={fmt(stats.cd, 2)} hint="per unit depth" />
          <Stat label="Cl" value={fmt(stats.cl, 2)} hint="per unit depth" />
        </dl>
      </div>

      {stats.warning ? (
        <p className="rounded-lg border border-border bg-surface px-3 py-2 text-sm text-accent">{stats.warning}</p>
      ) : null}

      <div className="flex flex-wrap gap-2">
        {MODES.map((item) => (
          <button
            key={item.id}
            type="button"
            aria-pressed={mode === item.id}
            onClick={() => setMode(item.id)}
            className={
              "min-h-11 rounded-lg border px-3 text-sm " +
              (mode === item.id
                ? "border-primary bg-surface-2 text-fg"
                : "border-border bg-surface text-muted")
            }
          >
            {item.label}
          </button>
        ))}
      </div>

      <div
        className={
          "h-2 rounded-lg " +
          (legendKind === "pressure" ? "legend-pressure" : legendKind === "viridis" ? "legend-viridis" : "legend-contrast")
        }
        aria-hidden="true"
      />
      <p className="font-mono text-xs text-muted">
        {legendLabel}
        {mode === "pressure"
          ? ` · ${fmt(stats.pressureMin, 4)} to ${fmt(stats.pressureMax, 4)} Pa`
          : mode === "vorticity"
            ? ` · ±${fmt(stats.vorticityMax, 2)} 1/s`
            : ` · 0 to ${fmt(stats.maxSpeed, 3)} m/s`}
        {` · t ${fmt(stats.time, 2)} s · div ${fmt(stats.divergence, 2)} · residual ${fmt(stats.residual, 2)}`}
      </p>

      <section className="grid gap-3 rounded-lg border border-border bg-surface p-3 sm:grid-cols-2">
        <label className="flex flex-col gap-1 text-sm text-muted">
          <span className="flex justify-between">
            Inlet speed <span className="font-mono text-fg">{fmt(spec.inlet, 3)} m/s</span>
          </span>
          <input
            className="lab-range"
            type="range"
            min={0.03}
            max={1.2}
            step={0.001}
            value={spec.inlet}
            onChange={(event) => applyInlet(Number(event.target.value))}
          />
        </label>
        <label className="flex flex-col gap-1 text-sm text-muted">
          <span className="flex justify-between">
            Angle <span className="font-mono text-fg">{fmt(spec.angleDeg, 0)}°</span>
          </span>
          <input
            className="lab-range"
            type="range"
            min={-25}
            max={35}
            step={1}
            value={spec.angleDeg}
            onChange={(event) => applyGeometry({ angleDeg: Number(event.target.value) })}
          />
        </label>
        <label className="flex flex-col gap-1 text-sm text-muted sm:col-span-2">
          <span className="flex justify-between">
            Size <span className="font-mono text-fg">{fmt(spec.diameter * 1000, 0)} mm</span>
          </span>
          <input
            className="lab-range"
            type="range"
            min={0.012}
            max={0.036}
            step={0.001}
            value={spec.diameter}
            onChange={(event) => applyGeometry({ diameter: Number(event.target.value) })}
          />
        </label>
      </section>

      <div className="flex flex-wrap gap-2">
        {BODIES.map((item) => (
          <button
            key={item.id}
            type="button"
            aria-pressed={spec.body === item.id}
            onClick={() => applyGeometry({ body: item.id })}
            className={
              "min-h-11 rounded-lg border px-3 text-sm " +
              (spec.body === item.id
                ? "border-primary bg-surface-2 text-fg"
                : "border-border bg-surface text-muted")
            }
          >
            {item.label}
          </button>
        ))}
        <button
          type="button"
          onClick={reset}
          className="inline-flex min-h-11 items-center gap-2 rounded-lg border border-border bg-surface px-3 text-sm text-fg"
        >
          <RotateCcw className="size-4" />
          Reset
        </button>
      </div>

      <div className="flex flex-wrap gap-2">
        {PRESETS.map((preset) => (
          <button
            key={preset.id}
            type="button"
            onClick={() =>
              applyGeometry({
                body: preset.body,
                inlet: preset.inlet,
                angleDeg: preset.angle,
                diameter: preset.diameter,
              })
            }
            className="min-h-11 rounded-lg border border-border bg-bg px-3 text-sm text-muted"
          >
            {preset.label}
          </button>
        ))}
      </div>

      <section className="rounded-lg border border-border bg-surface">
        <button
          type="button"
          className="flex min-h-11 w-full items-center justify-between px-3 text-left text-sm text-fg"
          aria-expanded={noteOpen}
          onClick={() => setNoteOpen((v) => !v)}
        >
          <span>What this solver is, and what it is not</span>
          <span className="font-mono text-xs text-muted">{noteOpen ? "Hide" : "Show"}</span>
        </button>
        {noteOpen ? (
          <div className="flex flex-col gap-3 border-t border-border px-3 py-3 text-sm text-pretty text-muted">
            <p>
              The native source is AeroLab 3D Beta (archive filename v0.8.0, notes still marked v0.7.0). It
              solves 3D incompressible Navier–Stokes on a staggered MAC grid with a PCG pressure solve,
              optional geometric multigrid, and an OpenGL particle view. Its own acceptance note says steady
              Cd, wake convergence, and the Mac frontend were not passed.
            </p>
            <p>
              This page keeps that tunnel setup in two dimensions: domain 0.20 m × 0.10 m, air ρ ={" "}
              {FLUID.density} kg/m³, μ = {FLUID.viscosity} Pa·s, inlet smoothstep over D/U, no-slip on the
              body, slip walls, and p = 0 on the outlet face. Advection is semi-Lagrangian so a phone-sized
              mesh can still grow a wake; the desktop code uses upwind or MUSCL. Pressure uses SOR, not the
              1e-7 PCG gate. Cd and Cl are section forces per unit depth from pressure and a wall shear
              estimate. They are not a certified drag coefficient.
            </p>
            <p>
              The nominal cylinder case is the AeroLab point: D = 20 mm, U = 0.113 m/s, Re ≈ 150. Give the
              ramp a second. Vortex shedding, if the mesh holds it, shows up in vorticity and as a wandering
              Cl.
            </p>
          </div>
        ) : null}
      </section>
    </div>
  );
}

function Stat({ label, value, hint }: { label: string; value: string; hint: string }) {
  return (
    <div className="bg-surface px-3 py-2">
      <dt className="font-mono text-xs tracking-wide text-muted uppercase">{label}</dt>
      <dd className="font-mono text-lg text-fg tabular-nums">{value}</dd>
      <p className="text-xs text-muted">{hint}</p>
    </div>
  );
}

function paint(
  canvas: HTMLCanvasElement,
  tunnel: WindTunnel,
  mode: FieldMode,
  luts: { contrast: Uint8ClampedArray; viridis: Uint8ClampedArray; pressure: Uint8ClampedArray },
  fieldRef: { current: HTMLCanvasElement | null },
): void {
  const parent = canvas.parentElement;
  if (!parent) return;
  const rect = parent.getBoundingClientRect();
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const w = Math.max(2, Math.floor(rect.width * dpr));
  const h = Math.max(2, Math.floor(rect.height * dpr));
  if (canvas.width !== w || canvas.height !== h) {
    canvas.width = w;
    canvas.height = h;
  }
  const ctx = canvas.getContext("2d");
  if (!ctx) return;

  let field = fieldRef.current;
  if (!field) {
    field = document.createElement("canvas");
    fieldRef.current = field;
  }
  if (field.width !== tunnel.nx || field.height !== tunnel.ny) {
    field.width = tunnel.nx;
    field.height = tunnel.ny;
  }
  const fctx = field.getContext("2d");
  if (!fctx) return;
  const image = fctx.createImageData(tunnel.nx, tunnel.ny);
  const data = image.data;
  const solidInk = [12, 22, 31];

  if (mode === "smoke") {
    const scale = Math.max(tunnel.inlet * 1.35, 0.05);
    const lut = luts.contrast;
    for (let j = 0; j < tunnel.ny; j += 1) {
      for (let i = 0; i < tunnel.nx; i += 1) {
        const id = tunnel.ic(i, j);
        const row = tunnel.ny - 1 - j;
        const pix = (row * tunnel.nx + i) * 4;
        if (tunnel.solid[id]) {
          data[pix] = solidInk[0];
          data[pix + 1] = solidInk[1];
          data[pix + 2] = solidInk[2];
        } else {
          const u = 0.5 * (tunnel.U[tunnel.iu(i, j)] + tunnel.U[tunnel.iu(i + 1, j)]);
          const v = 0.5 * (tunnel.V[tunnel.iv(i, j)] + tunnel.V[tunnel.iv(i, j + 1)]);
          const idx = shade(Math.hypot(u, v), 0, scale, null);
          data[pix] = lut[idx * 3] * 0.35;
          data[pix + 1] = lut[idx * 3 + 1] * 0.35;
          data[pix + 2] = lut[idx * 3 + 2] * 0.35;
        }
        data[pix + 3] = 255;
      }
    }
  } else {
    const scalar = tunnel.fillScalar(mode);
    const lut = mode === "speed" ? luts.viridis : luts.pressure;
    for (let j = 0; j < tunnel.ny; j += 1) {
      for (let i = 0; i < tunnel.nx; i += 1) {
        const id = tunnel.ic(i, j);
        const row = tunnel.ny - 1 - j;
        const pix = (row * tunnel.nx + i) * 4;
        if (tunnel.solid[id]) {
          data[pix] = solidInk[0];
          data[pix + 1] = solidInk[1];
          data[pix + 2] = solidInk[2];
        } else {
          const idx = shade(scalar.values[id], scalar.min, scalar.max, scalar.reference);
          data[pix] = lut[idx * 3];
          data[pix + 1] = lut[idx * 3 + 1];
          data[pix + 2] = lut[idx * 3 + 2];
        }
        data[pix + 3] = 255;
      }
    }
  }
  fctx.putImageData(image, 0, 0);
  ctx.imageSmoothingEnabled = true;
  ctx.drawImage(field, 0, 0, w, h);

  if (mode === "smoke") {
    const lut = luts.contrast;
    const scale = Math.max(tunnel.inlet * 1.35, 0.05);
    const parts = tunnel.particles;
    const radius = Math.max(1.2, w / 280);
    for (let i = 0; i < parts.length; i += 1) {
      const p = parts[i];
      const idx = shade(p.u, 0, scale, null);
      ctx.fillStyle = `rgb(${lut[idx * 3]},${lut[idx * 3 + 1]},${lut[idx * 3 + 2]})`;
      const px = (p.x / tunnel.lx) * w;
      const py = h - (p.y / tunnel.ly) * h;
      ctx.fillRect(px, py, radius, radius);
    }
  }

  traceBody(ctx, tunnel, w, h);
  if (mode === "speed") drawStreamlines(ctx, tunnel, w, h);

  ctx.strokeStyle = "rgba(46,196,182,0.35)";
  ctx.lineWidth = Math.max(1, w / 500);
  ctx.strokeRect(0.5, 0.5, w - 1, h - 1);
}
