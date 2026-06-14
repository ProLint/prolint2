/** Time series heatmap visualizations. */

import { memo, useMemo } from 'react';
import { Box, Typography } from '@mui/material';
import FastHeatmap from '../../../components/FastHeatmap';
import { HeatmapSkeleton } from './Skeletons';
import { getUnitLabel } from '../utils';
import type { TimeSeriesData, ResidueTimeSeriesData, InteractionData } from '../types';

interface TimeSeriesHeatmapProps {
  data: TimeSeriesData | null;
  loading: boolean;
  interactionData: InteractionData | null;
  onRowClick: (rowLabel: string) => void;
}

export const TimeSeriesHeatmap = memo(function TimeSeriesHeatmap({
  data,
  loading,
  interactionData,
  onRowClick,
}: TimeSeriesHeatmapProps) {
  const { matrix, colLabels, sortedResidues, showAsTime, unitLabel } = useMemo(() => {
    if (!data || !data.contact_counts) {
      return { matrix: [], colLabels: [], sortedResidues: [], showAsTime: false, unitLabel: '' };
    }

    const { contact_counts, frames, query_residues } = data;

    if (Object.keys(contact_counts).length === 0 || frames.length === 0) {
      return { matrix: [], colLabels: [], sortedResidues: [], showAsTime: false, unitLabel: '' };
    }

    const normalizeBy = interactionData?.params?.normalize_by || 'counts';
    const units = interactionData?.params?.units || 'ns';
    const normFactor = interactionData?.params?.norm_factor || 1;
    const showAsTime = normalizeBy === 'actual_time';
    const unitLabel = getUnitLabel(units);

    const sortedResidues = [...query_residues].sort((a, b) => a - b);
    const matrix = sortedResidues.map((resid) => contact_counts[resid] || new Array(frames.length).fill(0));
    const colLabels = showAsTime ? frames.map((f: number) => (f * normFactor).toFixed(2)) : frames.map(String);

    return { matrix, colLabels, sortedResidues, showAsTime, unitLabel };
  }, [data, interactionData]);

  if (loading || !data || !data.contact_counts) {
    return <HeatmapSkeleton height={300} rows={10} />;
  }

  if (matrix.length === 0) {
    return <Typography color="text.secondary">No time series data available</Typography>;
  }

  return (
    <FastHeatmap
      data={matrix}
      rowLabels={sortedResidues.map(String)}
      colLabels={colLabels}
      colorScheme="viridis"
      height={Math.min(Math.max(sortedResidues.length * 20 + 80, 200), 400)}
      showRowLabels={sortedResidues.length <= 30}
      xAxisLabel={showAsTime ? `Time (${unitLabel})` : "Frame"}
      yAxisLabel="Residue"
      disableSampling={true}
      onRowClick={onRowClick}
    />
  );
});

interface DatabaseContactsHeatmapProps {
  data: ResidueTimeSeriesData | null;
  loading: boolean;
  interactionData: InteractionData | null;
  selectedDatabaseType: string | null;
  onRowClick: (rowLabel: string) => void;
}

export const DatabaseContactsHeatmap = memo(function DatabaseContactsHeatmap({
  data,
  loading,
  interactionData,
  selectedDatabaseType,
  onRowClick,
}: DatabaseContactsHeatmapProps) {
  const { matrix, colLabels, sortedDbIds, showAsTime, unitLabel, rowTotals } = useMemo(() => {
    if (!data || !data.contact_matrix) {
      return { matrix: [], colLabels: [], sortedDbIds: [], showAsTime: false, unitLabel: '', rowTotals: [] };
    }

    const { contact_matrix, frames, database_ids } = data;

    if (database_ids.length === 0 || frames.length === 0) {
      return { matrix: [], colLabels: [], sortedDbIds: [], showAsTime: false, unitLabel: '', rowTotals: [] };
    }

    const normalizeBy = interactionData?.params?.normalize_by || 'counts';
    const units = interactionData?.params?.units || 'ns';
    const normFactor = interactionData?.params?.norm_factor || 1;
    const showAsTime = normalizeBy === 'actual_time';
    const unitLabel = getUnitLabel(units);

    const sortedDbIds = [...database_ids].sort((a, b) => a - b);
    const matrix = sortedDbIds.map((dbId) => contact_matrix[dbId] || new Array(frames.length).fill(0));
    const colLabels = showAsTime ? frames.map((f: number) => (f * normFactor).toFixed(2)) : frames.map(String);
    const rowTotals = matrix.map((row) => row.reduce((sum, val) => sum + val, 0));

    return { matrix, colLabels, sortedDbIds, showAsTime, unitLabel, rowTotals };
  }, [data, interactionData]);

  if (loading || !data || !data.contact_matrix) {
    return <HeatmapSkeleton height={300} rows={8} />;
  }

  if (matrix.length === 0) {
    return (
      <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: 200 }}>
        <Typography color="text.secondary">No contacts found for this residue.</Typography>
      </Box>
    );
  }

  const totalCount = data.total_database_ids;
  const displayedCount = sortedDbIds.length;
  const isTruncated = totalCount > displayedCount;

  return (
    <Box>
      <FastHeatmap
        data={matrix}
        rowLabels={sortedDbIds.map(String)}
        colLabels={colLabels}
        colorScheme="blues"
        height={Math.min(Math.max(sortedDbIds.length * 18 + 80, 200), 400)}
        showRowLabels={sortedDbIds.length <= 40}
        showColorbar={false}
        xAxisLabel={showAsTime ? `Time (${unitLabel})` : "Frame"}
        yAxisLabel={`${selectedDatabaseType} ID`}
        disableSampling={true}
        onRowClick={onRowClick}
        rowTotals={rowTotals}
      />
      {isTruncated && (
        <Typography variant="caption" color="text.secondary" sx={{ mt: 1, display: 'block' }}>
          Showing top {displayedCount} of the {totalCount} {selectedDatabaseType} molecules that had at least one contact.
        </Typography>
      )}
    </Box>
  );
});
