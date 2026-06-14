/** Shared contacts matrix and network visualizations. */

import { memo, useMemo } from 'react';
import { Box, Typography } from '@mui/material';
import ReactECharts from 'echarts-for-react';
import SharedContactsNetwork from '../../../components/SharedContactsNetwork';
import { colorScales } from '../../../theme/visualizationTheme';
import { HeatmapSkeleton } from './Skeletons';
import type { SharedContactsData, NetworkNode, NetworkEdge } from '../types';

interface SharedContactsMatrixProps {
  data: SharedContactsData | null;
  loading: boolean;
}

export const SharedContactsMatrix = memo(function SharedContactsMatrix({
  data,
  loading,
}: SharedContactsMatrixProps) {
  const option = useMemo(() => {
    if (!data || !data.matrix || data.matrix.length === 0) {
      return null;
    }

    const { labels, matrix } = data;
    const maxVal = Math.max(...matrix.flat(), 1);

    const echartsData: [number, number, number][] = [];
    matrix.forEach((row, yIdx) => {
      row.forEach((value, xIdx) => {
        echartsData.push([xIdx, yIdx, value]);
      });
    });

    return {
      tooltip: {
        trigger: 'item',
        formatter: (params: any) => `Res ${labels[params.data[0]]} ↔ Res ${labels[params.data[1]]}<br/>Shared: ${params.data[2]}`,
      },
      grid: { left: 10, right: 5, top: 5, bottom: 50, containLabel: false },
      xAxis: {
        type: 'category',
        data: labels.map(String),
        name: 'Query Residue',
        nameLocation: 'center',
        nameGap: 25,
        nameTextStyle: { fontSize: 11, fontWeight: 500 },
        axisLabel: { interval: Math.max(0, Math.floor(labels.length / 10) - 1), fontSize: 9, rotate: 0 },
        axisLine: { show: false, onZero: false },
        splitArea: { show: false },
      },
      yAxis: {
        type: 'category',
        data: labels.map(String),
        name: 'Query Residue',
        nameLocation: 'center',
        nameGap: 20,
        nameTextStyle: { fontSize: 11, fontWeight: 500 },
        axisLabel: { interval: Math.max(0, Math.floor(labels.length / 10) - 1), fontSize: 9 },
        splitArea: { show: false },
      },
      visualMap: {
        min: 0,
        max: maxVal,
        calculable: true,
        orient: 'horizontal',
        left: 'center',
        bottom: -16,
        itemHeight: 250,
        itemWidth: 7,
        textStyle: { fontSize: 9 },
        inRange: { color: colorScales.blues },
      },
      series: [{
        type: 'heatmap',
        data: echartsData,
        progressive: 0,
        progressiveThreshold: Infinity,
        animation: false,
      }],
    };
  }, [data]);

  if (loading) return <HeatmapSkeleton height={300} rows={8} />;

  if (!data || data.labels.length < 2) {
    return (
      <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '100%' }}>
        <Typography color="text.secondary">Not enough query residues with contacts to show matrix.</Typography>
      </Box>
    );
  }

  if (!option) return <Typography color="text.secondary">No data available</Typography>;

  return (
    <Box sx={{ width: '100%', height: '100%' }}>
      <ReactECharts option={option} style={{ width: '100%', height: '100%' }} opts={{ renderer: 'canvas' }} notMerge={true} />
    </Box>
  );
});

interface SharedContactsNetworkViewProps {
  data: SharedContactsData | null;
  loading: boolean;
  networkData: { nodes: NetworkNode[]; edges: NetworkEdge[] };
  selectedResiduesCount: number;
}

export const SharedContactsNetworkView = memo(function SharedContactsNetworkView({
  data,
  loading,
  networkData,
  selectedResiduesCount,
}: SharedContactsNetworkViewProps) {
  if (loading) return <HeatmapSkeleton height={300} rows={6} />;

  if (!data || data.labels.length < 2) {
    return (
      <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '100%', width: '100%' }}>
        <Typography color="text.secondary">Not enough query residues with contacts to show network.</Typography>
      </Box>
    );
  }

  if (selectedResiduesCount === 0) {
    return null;
  }

  return (
    <Box sx={{ width: '100%', height: '100%', minHeight: { xs: 220, sm: 300 } }}>
      <SharedContactsNetwork nodes={networkData.nodes} edges={networkData.edges} />
    </Box>
  );
});
