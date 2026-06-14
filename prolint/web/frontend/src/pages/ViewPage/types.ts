/** Type definitions for ViewPage. */

export interface InteractionData {
  composition: {
    resname_counts: Record<string, number>;
  };
  universe: {
    n_frames: number;
  };
  frame_range: {
    start: number;
    end: number;
    step: number;
  };
  params?: {
    units: string;
    normalize_by: string;
    norm_factor: number;
  };
}

export interface DensityMapData {
  x_edges: number[];
  y_edges: number[];
  density: number[][];
  query_density: number[][];
}

export interface LogoPlotData {
  residues: { resid: number; resname: string; chainID: string; value: number }[];
}

export interface SharedContactsData {
  labels: number[];
  matrix: number[][];
}

export interface TimeSeriesData {
  query_residues: number[];
  frames: number[];
  contact_counts: Record<string, number[]>;
}

export interface ResidueTimeSeriesData {
  database_ids: number[];
  frames: number[];
  contact_matrix: Record<string, number[]>;
  total_database_ids: number;
}

export interface DistanceTimeData {
  frames: number[];
  distances: number[];
  min_distances: number[];
  contact_frames: number[];
  positions: {
    query: { x: number; y: number }[];
    database: { x: number; y: number }[];
  };
}

export interface AtomDistancesData {
  frame: number;
  query_atoms: string[];
  database_atoms: string[];
  distance_matrix: number[][];
  min_distance: number;
  max_distance: number;
}

export interface MonoFit {
  k_off: number;
  r_squared: number;
  aic: number;
  fitted_curve: number[];
  half_life: number | null;
}

export interface BiFit {
  a_fast: number;
  k_fast: number;
  k_slow: number;
  r_squared: number;
  aic: number;
  fitted_curve: number[];
  half_life_fast: number | null;
  half_life_slow: number | null;
}

export interface KineticsData {
  mode: 'individual' | 'accumulated';
  kinetics: {
    koff: number;
    kon: number;
    kd: number | null;
    mean_residence_time: number;
    std_residence_time: number;
    max_residence_time: number;
    occupancy: number;
    n_events: number;
    total_contact_frames: number;
    total_trajectory_frames: number;
  };
  survival_curve: {
    lag_times: number[];
    survival_probability: number[];
    mono_fit: MonoFit | null;
    bi_fit: BiFit | null;
    selected_model: 'monoexponential' | 'biexponential' | null;
    min_events_mono: number;
    min_events_bi: number;
  };
  residence_distribution: {
    bins: number[];
    counts: number[];
  };
  events: { start: number; end: number; duration: number }[];
}

export interface CompositionItem {
  name: string;
  value: number;
  percentage: number;
  color: string;
}

export interface ExpandedSections {
  overview: boolean;
  residueSelection: boolean;
  timeSeries: boolean;
  pairAnalysis: boolean;
}

export interface NetworkNode {
  id: string;
  label: string;
  type: string;
  residue_id: number;
  restype: string;
}

export interface NetworkEdge {
  source: string;
  target: string;
  value: number;
}
