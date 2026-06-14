/** Shared color utilities. Loads theme from theme.json (single source of truth). */

// Import the shared theme JSON
// Note: This path is resolved via Vite alias '@prolint/config'
import sharedTheme from '@prolint/config/theme.json';

// =============================================================================
// Type Definitions
// =============================================================================

export interface RGBA {
  r: number;
  g: number;
  b: number;
  a: number;
}

export interface GradientStop {
  position: number;
  color: RGBA;
}

export interface ThemeColors {
  primary: Record<string, string>;
  accent: Record<string, string>;
  success: Record<string, string>;
  warning: Record<string, string>;
  error: Record<string, string>;
  info: Record<string, string>;
  neutral: Record<string, string>;
  data: Record<string, string>;
  background: Record<string, string>;
  text: Record<string, string>;
  border: Record<string, string>;
}

// =============================================================================
// Theme Exports (from shared JSON)
// =============================================================================

export const COLORS = sharedTheme.colors as ThemeColors;
export const COLOR_SCALES = sharedTheme.colorScales as Record<string, string[]>;
export const GRADIENTS = sharedTheme.gradients as Record<string, GradientStop[]>;
export const AMINO_ACID_COLORS = sharedTheme.aminoAcidColors as Record<string, string>;
export const AMINO_ACID_ONE_LETTER = sharedTheme.aminoAcidOneLetter as Record<string, string>;
export const UNIT_LABELS = sharedTheme.unitLabels as Record<string, string>;

// =============================================================================
// Color Conversion Utilities
// =============================================================================

/**
 * Convert hex color to RGB tuple.
 *
 * @param hexColor - Hex color string like '#06b6d4' or '06b6d4'
 * @returns Object with r, g, b values (0-255)
 *
 * @example
 * hexToRgb('#06b6d4') // { r: 6, g: 182, b: 212 }
 */
export function hexToRgb(hexColor: string): { r: number; g: number; b: number } {
  const hex = hexColor.replace('#', '');
  return {
    r: parseInt(hex.substring(0, 2), 16),
    g: parseInt(hex.substring(2, 4), 16),
    b: parseInt(hex.substring(4, 6), 16),
  };
}

/**
 * Convert RGB values to hex color.
 *
 * @param r - Red value (0-255)
 * @param g - Green value (0-255)
 * @param b - Blue value (0-255)
 * @returns Hex color string like '#06b6d4'
 *
 * @example
 * rgbToHex(6, 182, 212) // '#06b6d4'
 */
export function rgbToHex(r: number, g: number, b: number): string {
  const toHex = (n: number) => n.toString(16).padStart(2, '0');
  return `#${toHex(r)}${toHex(g)}${toHex(b)}`;
}

/**
 * Convert RGBA values to CSS rgba() string.
 *
 * @param r - Red value (0-255)
 * @param g - Green value (0-255)
 * @param b - Blue value (0-255)
 * @param a - Alpha value (0-1)
 * @returns CSS rgba string like 'rgba(6, 182, 212, 0.7)'
 *
 * @example
 * rgbaToCss(6, 182, 212, 0.7) // 'rgba(6, 182, 212, 0.7)'
 */
export function rgbaToCss(r: number, g: number, b: number, a: number): string {
  return `rgba(${r}, ${g}, ${b}, ${a})`;
}

/**
 * Parse CSS rgba() string to object.
 *
 * @param rgbaStr - CSS rgba string like 'rgba(6, 182, 212, 0.7)'
 * @returns Object with r, g, b (0-255) and a (0-1)
 *
 * @example
 * cssRgbaToObject('rgba(6, 182, 212, 0.7)') // { r: 6, g: 182, b: 212, a: 0.7 }
 */
export function cssRgbaToObject(rgbaStr: string): RGBA {
  const match = rgbaStr.match(/rgba?\((\d+),\s*(\d+),\s*(\d+),?\s*([\d.]+)?\)/);
  if (match) {
    return {
      r: parseInt(match[1], 10),
      g: parseInt(match[2], 10),
      b: parseInt(match[3], 10),
      a: match[4] ? parseFloat(match[4]) : 1.0,
    };
  }
  return { r: 128, g: 128, b: 128, a: 1.0 };
}

// =============================================================================
// Gradient Interpolation
// =============================================================================

/**
 * Get a specific color stop from a gradient.
 *
 * @param gradientName - Name of the gradient: 'sharedContacts' or 'density'
 * @param position - Index of the color stop (0, 1, 2, ...)
 * @returns CSS rgba string
 */
export function getGradientColor(gradientName: string, position: number): string {
  const gradient = GRADIENTS[gradientName] || GRADIENTS.sharedContacts;
  const idx = Math.min(position, gradient.length - 1);
  const stop = gradient[idx];
  return rgbaToCss(stop.color.r, stop.color.g, stop.color.b, stop.color.a);
}

/**
 * Interpolate a color from a gradient based on a value.
 *
 * @param value - The value to map to a color
 * @param minVal - Minimum value in range
 * @param maxVal - Maximum value in range
 * @param gradientName - Name of the gradient: 'sharedContacts' or 'density'
 * @returns CSS rgba string
 *
 * @example
 * interpolateGradient(50, 0, 100, 'sharedContacts') // 'rgba(168, 85, 247, 0.7)'
 */
export function interpolateGradient(
  value: number,
  minVal: number,
  maxVal: number,
  gradientName: string = 'sharedContacts'
): string {
  const gradient = GRADIENTS[gradientName] || GRADIENTS.sharedContacts;

  // Normalize value to 0-1 range
  let normalized = maxVal > minVal ? (value - minVal) / (maxVal - minVal) : 0.5;
  normalized = Math.max(0, Math.min(1, normalized));

  // Find the two stops to interpolate between
  for (let i = 0; i < gradient.length - 1; i++) {
    const pos1 = gradient[i].position;
    const pos2 = gradient[i + 1].position;

    if (pos1 <= normalized && normalized <= pos2) {
      const t = pos2 > pos1 ? (normalized - pos1) / (pos2 - pos1) : 0;

      const c1 = gradient[i].color;
      const c2 = gradient[i + 1].color;

      const r = Math.round(c1.r + (c2.r - c1.r) * t);
      const g = Math.round(c1.g + (c2.g - c1.g) * t);
      const b = Math.round(c1.b + (c2.b - c1.b) * t);
      const a = c1.a + (c2.a - c1.a) * t;

      return rgbaToCss(r, g, b, a);
    }
  }

  // Fallback to last color
  const last = gradient[gradient.length - 1].color;
  return rgbaToCss(last.r, last.g, last.b, last.a);
}

/**
 * Get a hex color for a value from a discrete color scale.
 *
 * @param value - The value to map
 * @param scaleName - Name of the color scale
 * @param vmin - Minimum value
 * @param vmax - Maximum value
 * @returns Hex color string
 *
 * @example
 * getColorForValue(0.5, 'prolint', 0, 1) // '#22d3ee'
 */
export function getColorForValue(
  value: number,
  scaleName: string = 'prolint',
  vmin: number = 0,
  vmax: number = 1
): string {
  const colors = COLOR_SCALES[scaleName] || COLOR_SCALES.prolint;
  let normalized = vmax > vmin ? (value - vmin) / (vmax - vmin) : 0.5;
  normalized = Math.max(0, Math.min(1, normalized));
  const idx = Math.floor(normalized * (colors.length - 1));
  return colors[idx];
}

/**
 * Get display label for a time unit.
 *
 * @param unit - Time unit: 'fs', 'ps', 'ns', 'us', 'ms', 's'
 * @returns Display label (e.g., 'us' -> 'μs')
 */
export function getUnitLabel(unit: string): string {
  return UNIT_LABELS[unit] || unit;
}

// =============================================================================
// Legacy gradient format (for backwards compatibility)
// =============================================================================

/**
 * Convert JSON gradient to legacy [position, cssString][] format.
 */
function convertGradientToLegacy(gradient: GradientStop[]): [number, string][] {
  return gradient.map((stop) => [
    stop.position,
    rgbaToCss(stop.color.r, stop.color.g, stop.color.b, stop.color.a),
  ]);
}

export const SHARED_CONTACTS_GRADIENT = convertGradientToLegacy(GRADIENTS.sharedContacts);
export const DENSITY_GRADIENT = convertGradientToLegacy(GRADIENTS.density);
