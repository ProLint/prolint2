/**
 * Type declarations for @prolint/config module
 */

declare module '@prolint/config/theme.json' {
  interface RGBA {
    r: number;
    g: number;
    b: number;
    a: number;
  }

  interface GradientStop {
    position: number;
    color: RGBA;
  }

  interface Theme {
    _comment: string;
    colors: {
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
    };
    colorScales: {
      mako: string[];
      viridis: string[];
      prolint: string[];
      blues: string[];
      diverging: string[];
      categorical: string[];
    };
    gradients: {
      sharedContacts: GradientStop[];
      density: GradientStop[];
    };
    aminoAcidColors: Record<string, string>;
    aminoAcidOneLetter: Record<string, string>;
    unitLabels: Record<string, string>;
  }

  const theme: Theme;
  export default theme;
}
