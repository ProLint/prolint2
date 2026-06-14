/** Color utility functions for visualizations. */

/** Interpolate between two RGBA colors. */
export function interpolateRgba(color1: string, color2: string, t: number): string {
  const match1 = color1.match(/rgba?\((\d+),\s*(\d+),\s*(\d+),?\s*([\d.]+)?\)/);
  const match2 = color2.match(/rgba?\((\d+),\s*(\d+),\s*(\d+),?\s*([\d.]+)?\)/);

  if (!match1 || !match2) return color1;

  const r = Math.round(parseInt(match1[1]) + (parseInt(match2[1]) - parseInt(match1[1])) * t);
  const g = Math.round(parseInt(match1[2]) + (parseInt(match2[2]) - parseInt(match1[2])) * t);
  const b = Math.round(parseInt(match1[3]) + (parseInt(match2[3]) - parseInt(match1[3])) * t);
  const a = parseFloat(match1[4] || '1');

  return `rgba(${r}, ${g}, ${b}, ${a})`;
}

/** Interpolate a color along a multi-point color scale. */
export function interpolateColorScale(
  value: number,
  min: number,
  max: number,
  scale: [number, string][]
): string {
  const normalized = max > min ? (value - min) / (max - min) : 0.5;

  if (normalized <= 0) return String(scale[0][1]);
  if (normalized >= 1) return String(scale[scale.length - 1][1]);

  // Find the two colors to interpolate between
  for (let i = 0; i < scale.length - 1; i++) {
    const [pos1, color1] = scale[i];
    const [pos2, color2] = scale[i + 1];

    if (normalized >= pos1 && normalized <= pos2) {
      const t = (normalized - pos1) / (pos2 - pos1);
      return interpolateRgba(String(color1), String(color2), t);
    }
  }

  return String(scale[scale.length - 1][1]);
}

/** Convert hex color to RGB object. */
export function hexToRgb(hex: string): { r: number; g: number; b: number } | null {
  const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  return result
    ? {
        r: parseInt(result[1], 16),
        g: parseInt(result[2], 16),
        b: parseInt(result[3], 16),
      }
    : null;
}

/** Calculate relative luminance for contrast checking. */
export function getRelativeLuminance(r: number, g: number, b: number): number {
  const [rs, gs, bs] = [r, g, b].map((c) => {
    const s = c / 255;
    return s <= 0.03928 ? s / 12.92 : Math.pow((s + 0.055) / 1.055, 2.4);
  });
  return 0.2126 * rs + 0.7152 * gs + 0.0722 * bs;
}

/** Determine if text should be light or dark based on background (WCAG). */
export function getContrastTextColor(backgroundColor: string): 'light' | 'dark' {
  const rgb = hexToRgb(backgroundColor);
  if (!rgb) return 'dark';

  const luminance = getRelativeLuminance(rgb.r, rgb.g, rgb.b);
  // WCAG recommends 4.5:1 contrast ratio for normal text
  // Luminance > 0.179 means background is light, so use dark text
  return luminance > 0.179 ? 'dark' : 'light';
}
