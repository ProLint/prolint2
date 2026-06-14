/** ProLint visualization theme. Colors loaded from shared theme.json. */

import {
  COLORS,
  COLOR_SCALES,
  AMINO_ACID_COLORS,
  SHARED_CONTACTS_GRADIENT,
  DENSITY_GRADIENT,
  hexToRgb,
  rgbaToCss,
} from './sharedColors';

// ============================================================================
// Legacy colors export (for backwards compatibility)
// ============================================================================

export const colors = {
  primary: COLORS.primary,
  accent: COLORS.accent,
  success: COLORS.success,
  warning: COLORS.warning,
  error: COLORS.error,
  info: COLORS.info,
  neutral: COLORS.neutral,
  data: {
    query: COLORS.data.query,
    queryDark: COLORS.data.query_dark,
    queryLight: COLORS.data.query_light,
    database: COLORS.data.database,
    databaseDark: COLORS.data.database_dark,
    databaseLight: COLORS.data.database_light,
    highlight: COLORS.data.highlight,
    highlightDark: COLORS.data.highlight_dark,
  },
  background: COLORS.background,
  text: COLORS.text,
  border: COLORS.border,
} as const;

// ============================================================================
// Scientific Color Scales for Data Visualization
// ============================================================================

export const colorScales = {
  // Sequential scales from shared
  mako: COLOR_SCALES.mako,
  viridis: COLOR_SCALES.viridis,
  prolint: COLOR_SCALES.prolint,
  blues: COLOR_SCALES.blues,
  diverging: COLOR_SCALES.diverging,
  categorical: COLOR_SCALES.categorical,

  // Shared contacts gradient (connection strength)
  sharedContacts: SHARED_CONTACTS_GRADIENT,

  // Density gradient
  density: DENSITY_GRADIENT,
} as const;

// ============================================================================
// Visualization Theme Export
// ============================================================================

export const visualizationTheme = {
  // Flattened colors for backwards compatibility
  colors: {
    // Primary brand colors
    primary: colors.primary['500'],
    primaryDark: colors.primary['700'],
    primaryLight: colors.primary['400'],

    // Accent colors
    accent: colors.accent['500'],
    accentDark: colors.accent['600'],
    accentLight: colors.accent['400'],

    // Semantic colors
    success: colors.success.main,
    warning: colors.warning.main,
    error: colors.error.main,
    info: colors.info.main,

    // Network/graph colors
    nodeQuery: colors.data.query,
    nodeQueryFocus: colors.data.highlight,
    nodeDatabase: colors.data.database,
    edgeDefault: (() => {
      const { r, g, b } = hexToRgb(colors.neutral['500']);
      return rgbaToCss(r, g, b, 0.3);
    })(),

    // Heatmap color scales
    heatmapSequential: [colors.primary['100'], colors.primary['800']],
    heatmapDiverging: [colors.primary['700'], colors.neutral['100'], colors.error.dark],

    // Background
    background: colors.background.default,
    backgroundAlt: colors.background.subtle,
    surface: colors.background.paper,

    // Text
    textPrimary: colors.text.primary,
    textSecondary: colors.text.secondary,
    textTertiary: colors.text.tertiary,

    // Borders
    border: colors.border.default,
    borderLight: colors.border.light,
  },

  // Color scales for Plotly and other visualizations
  colorscales: {
    viridis: 'Viridis',
    plasma: 'Plasma',
    inferno: 'Inferno',
    turbo: 'Turbo',
    blues: 'Blues',
    rdbu: 'RdBu',
    sharedContacts: colorScales.sharedContacts,
    density: colorScales.density,
    categorical: colorScales.categorical,
    prolint: colorScales.prolint,
    mako: colorScales.mako,
  },

  // Typography
  typography: {
    fontFamily: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif',
    fontSize: {
      small: 10,
      normal: 12,
      medium: 14,
      large: 16,
      xlarge: 18,
    },
    fontWeight: {
      normal: 400,
      medium: 500,
      semibold: 600,
      bold: 700,
    },
  },

  // Layout spacing
  spacing: {
    xs: 4,
    sm: 8,
    md: 16,
    lg: 24,
    xl: 32,
    xxl: 48,
  },

  // Plot dimensions
  plotDimensions: {
    barChartHeight: 400,
    heatmapMinHeight: 300,
    heatmapMaxHeight: 800,
    networkHeight: 600,
    densityMapHeight: 500,
  },

  // Common Plotly layout config
  plotlyLayout: {
    font: {
      family: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif',
      size: 12,
      color: colors.text.primary,
    },
    paper_bgcolor: colors.background.default,
    plot_bgcolor: colors.background.subtle,
    hovermode: 'closest' as const,
    hoverlabel: {
      bgcolor: colors.background.default,
      bordercolor: colors.border.default,
      font: {
        family: '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif',
        size: 11,
        color: colors.text.primary,
      },
    },
    margin: {
      l: 60,
      r: 120,
      t: 50,
      b: 60,
    },
  },

  // Common Plotly config
  plotlyConfig: {
    responsive: true,
    displaylogo: false,
    displayModeBar: true,
    modeBarButtonsToRemove: ['lasso2d', 'select2d'] as string[],
    toImageButtonOptions: {
      format: 'png' as const,
      filename: 'prolint_plot',
      height: 1200,
      width: 1600,
      scale: 2,
    },
  },

  // Network (Cytoscape) styles
  cytoscapeStyle: [
    {
      selector: 'node',
      style: {
        'background-color': colors.data.query,
        'label': 'data(label)',
        'color': colors.text.primary,
        'text-valign': 'center',
        'text-halign': 'center',
        'font-size': '10px',
        'font-family': '-apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif',
        'font-weight': '600',
        'text-outline-width': 2,
        'text-outline-color': colors.background.default,
        'width': 30,
        'height': 30,
        'border-width': 2,
        'border-color': colors.data.queryDark,
      },
    },
    {
      selector: 'node:selected',
      style: {
        'background-color': colors.data.highlight,
        'border-color': colors.data.highlightDark,
        'border-width': 3,
        'width': 40,
        'height': 40,
      },
    },
    {
      selector: 'node.focused',
      style: {
        'background-color': colors.data.highlight,
        'border-color': colors.data.highlightDark,
        'border-width': 3,
        'width': 40,
        'height': 40,
      },
    },
    {
      selector: 'node.neighbor',
      style: {
        'background-color': colors.data.queryLight,
        'border-color': colors.data.query,
      },
    },
    {
      selector: 'edge',
      style: {
        'width': 'data(width)',
        'line-color': 'data(color)',
        'curve-style': 'bezier',
        'opacity': 0.7,
      },
    },
    {
      selector: 'edge:selected',
      style: {
        'opacity': 1,
        'width': 4,
      },
    },
  ],

  // Amino acid colors (from shared)
  aminoAcidColors: AMINO_ACID_COLORS,

  // Transitions
  transitions: {
    fast: '150ms cubic-bezier(0.4, 0, 0.2, 1)',
    normal: '300ms cubic-bezier(0.4, 0, 0.2, 1)',
    slow: '500ms cubic-bezier(0.4, 0, 0.2, 1)',
  },
};

export type VisualizationTheme = typeof visualizationTheme;
export type Colors = typeof colors;
export type ColorScales = typeof colorScales;
