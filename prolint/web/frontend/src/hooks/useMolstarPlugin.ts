/** Shared hook for Molstar plugin initialization with minimal config for fast loading. */

import { useEffect, useRef, useState, useCallback } from 'react';
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import { renderReact18 } from 'molstar/lib/mol-plugin-ui/react18';
import { PluginUISpec } from 'molstar/lib/mol-plugin-ui/spec';
import { PluginUIContext } from 'molstar/lib/mol-plugin-ui/context';
import { PluginConfig } from 'molstar/lib/mol-plugin/config';
import { PluginSpec } from 'molstar/lib/mol-plugin/spec';
import { PluginBehaviors } from 'molstar/lib/mol-plugin/behavior';
import { StructureElement } from 'molstar/lib/mol-model/structure';
import { Loci } from 'molstar/lib/mol-model/loci';
import { OrderedSet } from 'molstar/lib/mol-data/int';
import type { LociLabelProvider } from 'molstar/lib/mol-plugin-state/manager/loci-label';

export interface UseMolstarPluginResult {
  containerRef: React.RefObject<HTMLDivElement>;
  plugin: PluginUIContext | null;
  isReady: boolean;
  isLoading: boolean;
  error: string | null;
  loadIdRef: React.MutableRefObject<number>;
  setIsLoading: (loading: boolean) => void;
  setError: (error: string | null) => void;
  loadPdbData: (pdbData: string) => Promise<{
    data: unknown;
    trajectory: unknown;
    model: unknown;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    structure: any;
  } | null>;
}

// Minimal plugin spec for fast initialization with hover support
const MinimalPluginSpec: PluginUISpec = {
  actions: [],
  behaviors: [
    // Essential behaviors for hover/highlight functionality
    PluginSpec.Behavior(PluginBehaviors.Representation.HighlightLoci),
    // Note: We add a custom label provider after init instead of DefaultLociLabelProvider
    PluginSpec.Behavior(PluginBehaviors.Camera.FocusLoci),
    PluginSpec.Behavior(PluginBehaviors.Camera.CameraControls),
  ],
  animations: [],
  customParamEditors: [],
  layout: {
    initial: {
      isExpanded: false,
      showControls: false,
      controlsDisplay: 'reactive',
    },
  },
  components: {
    remoteState: 'none',
  },
  config: [
    [PluginConfig.VolumeStreaming.Enabled, false],
    [PluginConfig.Viewport.ShowExpand, false],
    [PluginConfig.Viewport.ShowControls, false],
    [PluginConfig.Viewport.ShowSettings, false],
    [PluginConfig.Viewport.ShowSelectionMode, false],
    [PluginConfig.Viewport.ShowAnimation, false],
  ],
};

/** Label provider showing residue type and ID on hover. */
const residueLabelProvider: LociLabelProvider = {
  label: (loci: Loci): string | undefined => {
    if (!StructureElement.Loci.is(loci)) return undefined;
    if (StructureElement.Loci.isEmpty(loci)) return undefined;

    const { unit, indices } = loci.elements[0];
    if (OrderedSet.size(indices) === 0) return undefined;

    const elementIndex = OrderedSet.getAt(indices, 0);
    const atomIndex = unit.elements[elementIndex];
    const residueIndex = unit.model.atomicHierarchy.residueAtomSegments.index[atomIndex];

    const residueName = unit.model.atomicHierarchy.atoms.label_comp_id.value(atomIndex);
    const residueId = unit.model.atomicHierarchy.residues.auth_seq_id.value(residueIndex);

    return `${residueName} ${residueId}`;
  },
  priority: 1,
};

export function useMolstarPlugin(): UseMolstarPluginResult {
  const containerRef = useRef<HTMLDivElement>(null);
  const pluginRef = useRef<PluginUIContext | null>(null);
  const loadIdRef = useRef(0);

  const [isReady, setIsReady] = useState(false);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Initialize Molstar plugin with minimal spec
  useEffect(() => {
    if (!containerRef.current || pluginRef.current) return;

    let cancelled = false;

    const initPlugin = async () => {
      try {
        const plugin = await createPluginUI({
          target: containerRef.current!,
          render: renderReact18,
          spec: MinimalPluginSpec,
        });

        if (cancelled) {
          plugin.dispose();
          return;
        }

        // Register custom label provider for residue type and ID only
        plugin.managers.lociLabels.addProvider(residueLabelProvider);

        pluginRef.current = plugin;
        setIsReady(true);
      } catch (err: unknown) {
        if (!cancelled) {
          const errorMessage = err instanceof Error ? err.message : String(err);
          setError('Failed to initialize 3D viewer: ' + errorMessage);
        }
      }
    };

    initPlugin();

    return () => {
      cancelled = true;
      if (pluginRef.current) {
        pluginRef.current.dispose();
        pluginRef.current = null;
      }
      setIsReady(false);
    };
  }, []);

  // Optimized PDB loading - minimal overhead
  const loadPdbData = useCallback(async (pdbData: string) => {
    const plugin = pluginRef.current;
    if (!plugin || !pdbData) return null;

    const currentLoadId = loadIdRef.current;

    // Clear and load in sequence
    await plugin.clear();
    if (loadIdRef.current !== currentLoadId) return null;

    const data = await plugin.builders.data.rawData({
      data: pdbData,
      label: 'structure',
    });

    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const model = await plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.builders.structure.createStructure(model);

    if (loadIdRef.current !== currentLoadId) return null;

    return { data, trajectory, model, structure };
  }, []);

  return {
    containerRef,
    plugin: pluginRef.current,
    isReady,
    isLoading,
    error,
    loadIdRef,
    setIsLoading,
    setError,
    loadPdbData,
  };
}

/** Convert hex color to Mol* Color integer. */
export function hexToMolColor(hex: string): number {
  return parseInt(hex.replace('#', ''), 16);
}
