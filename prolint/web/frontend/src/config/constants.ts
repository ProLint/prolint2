/** Application-wide constants and configuration values. */

export const LAYOUT = {
  DRAWER_WIDTH: 200,
  MIN_VIEWER_HEIGHT: 350,
} as const;

export const HEATMAP = {
  MAX_COLUMNS: 200,
  ROW_HEIGHT: {
    SMALL: 18,
    MEDIUM: 12,
    LARGE: 8,
  },
  HEIGHT_PADDING: 100,
  MAX_ROWS_FOR_LABELS: 40,
} as const;

export const NETWORK = {
  MAX_HEIGHT: 450,
  RESIZE_DELAY_MS: 200,
  ANIMATION_DURATION_MS: 300,
  MIN_ZOOM: 0.1,
  MAX_ZOOM: 3,
} as const;

export const UI = {
  CHIP: {
    HEIGHT_MOBILE: 28,
    HEIGHT_DESKTOP: 24,
    FONT_SIZE_MOBILE: '0.65rem',
    FONT_SIZE_DESKTOP: '0.7rem',
  },
} as const;

export const TIMEOUTS = {
  MOLSTAR_INIT_MS: 30000,
  STRUCTURE_LOAD_MS: 60000,
  RESIDUE_SELECTION_DEBOUNCE_MS: 300,
} as const;

export const FORMATTING = {
  EXPONENTIAL_PRECISION: 1,
  DECIMAL_PRECISION: 2,
  PERCENTAGE_PRECISION: 1,
} as const;
