/**
 * TypeScript type definitions for ProLint2 Dashboard
 */

// Replica Detection Types

export interface ReplicaInfo {
  replica_id: string;
  n_residues: number;
  first_resid: number;
  last_resid: number;
}

export interface ReplicaAnalysisResult {
  success: boolean;
  n_atoms: number;
  n_residues: number;
  n_replicas: number;
  detection_method?: 'bond_connectivity' | 'resid_replication' | null;
  has_repeated_resids?: boolean | null;
  replica_info?: ReplicaInfo[];
  message: string;
}

// Dataset Types

export interface DatasetInfo {
  id: string;
  name: string;
  topology_file?: string;
  trajectory_file?: string;
  n_frames?: number;
  n_atoms?: number;
  n_residues?: number;
  status: string;
}

// Computation Types

export interface ComputeResponse {
  result_id: string;
  status: string;
  computation_time: number;
}

// Interaction Data Types (returned from get_interaction_data)

export interface InteractionParams {
  units: string;
  normalize_by: string;
  norm_factor: number;
  selected_replica: string | null;
  replica_info: ReplicaInfo[];
}

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
  params: InteractionParams;
}
