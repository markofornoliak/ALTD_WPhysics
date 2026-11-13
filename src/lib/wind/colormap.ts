/** AeroLab 3D color stops (visualization/ColorMaps.cpp), kept in the same order. */

export type MapKind = "viridis" | "contrast" | "pressure" | "plasma";

type Stop = { t: number; r: number; g: number; b: number };

const PLASMA: Stop[] = [
  { t: 0, r: 0.05, g: 0.03, b: 0.528 },
  { t: 0.18, r: 0.365, g: 0.003, b: 0.65 },
  { t: 0.4, r: 0.692, g: 0.165, b: 0.565 },
  { t: 0.62, r: 0.902, g: 0.33, b: 0.329 },
  { t: 0.82, r: 0.991, g: 0.587, b: 0.12 },
  { t: 1, r: 0.94, g: 0.975, b: 0.131 },
];

const VIRIDIS: Stop[] = [
  { t: 0, r: 0.267, g: 0.005, b: 0.329 },
  { t: 0.25, r: 0.23, g: 0.322, b: 0.545 },
  { t: 0.5, r: 0.128, g: 0.567, b: 0.551 },
  { t: 0.75, r: 0.369, g: 0.789, b: 0.383 },
  { t: 1, r: 0.993, g: 0.906, b: 0.144 },
];

const CONTRAST: Stop[] = [
  { t: 0, r: 0.02, g: 0.055, b: 0.13 },
  { t: 0.22, r: 0, g: 0.42, b: 0.71 },
  { t: 0.48, r: 0, g: 0.85, b: 0.79 },
  { t: 0.72, r: 0.965, g: 0.78, b: 0.09 },
  { t: 1, r: 0.9, g: 0.08, b: 0.04 },
];

const PRESSURE: Stop[] = [
  { t: 0, r: 0.055, g: 0.12, b: 0.36 },
  { t: 0.22, r: 0.05, g: 0.43, b: 0.72 },
  { t: 0.43, r: 0.3, g: 0.79, b: 0.85 },
  { t: 0.5, r: 0.94, g: 0.945, b: 0.925 },
  { t: 0.63, r: 0.965, g: 0.65, b: 0.28 },
  { t: 0.82, r: 0.85, g: 0.22, b: 0.1 },
  { t: 1, r: 0.48, g: 0.025, b: 0.06 },
];

function stopsFor(kind: MapKind): Stop[] {
  if (kind === "viridis") return VIRIDIS;
  if (kind === "contrast") return CONTRAST;
  if (kind === "pressure") return PRESSURE;
  return PLASMA;
}

export function sampleColor(kind: MapKind, t: number): [number, number, number] {
  const stops = stopsFor(kind);
  const x = t < 0 ? 0 : t > 1 ? 1 : t;
  let hi = 1;
  while (hi < stops.length && x > stops[hi].t) hi += 1;
  if (hi >= stops.length) hi = stops.length - 1;
  const lo = hi > 0 ? hi - 1 : 0;
  const span = stops[hi].t - stops[lo].t;
  const f = span > 0 ? (x - stops[lo].t) / span : 0;
  const a = stops[lo];
  const b = stops[hi];
  return [a.r + (b.r - a.r) * f, a.g + (b.g - a.g) * f, a.b + (b.b - a.b) * f];
}

/** 256-entry RGB LUT, values 0–255. */
export function buildLut(kind: MapKind): Uint8ClampedArray {
  const lut = new Uint8ClampedArray(256 * 3);
  for (let i = 0; i < 256; i += 1) {
    const [r, g, b] = sampleColor(kind, i / 255);
    lut[i * 3] = Math.round(r * 255);
    lut[i * 3 + 1] = Math.round(g * 255);
    lut[i * 3 + 2] = Math.round(b * 255);
  }
  return lut;
}
