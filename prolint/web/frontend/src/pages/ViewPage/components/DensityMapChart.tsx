/** 2D density map visualization. */

import { memo, useMemo } from 'react';
import { Box, Alert, CircularProgress, Skeleton } from '@mui/material';
import ReactECharts from 'echarts-for-react';
import { colorScales } from '../../../theme/visualizationTheme';
import type { DensityMapData } from '../types';

interface DensityMapChartProps {
  densityMap: DensityMapData | null;
  loading: boolean;
  startFrame: number;
  endFrame: number;
}

const DensityMapChart = memo(function DensityMapChart({
  densityMap,
  loading,
  startFrame,
  endFrame,
}: DensityMapChartProps) {
  const option = useMemo(() => {
    if (!densityMap || !densityMap.density || densityMap.density.length === 0) {
      return null;
    }

    const xCenters = densityMap.x_edges.slice(0, -1).map((x: number, i: number) =>
      ((x + densityMap.x_edges[i + 1]) / 2).toFixed(1)
    );
    const yCenters = densityMap.y_edges.slice(0, -1).map((y: number, i: number) =>
      ((y + densityMap.y_edges[i + 1]) / 2).toFixed(1)
    );

    const echartsData: [number, number, number][] = [];
    let minVal = Infinity, maxVal = -Infinity;

    densityMap.density.forEach((row: number[], yIdx: number) => {
      row.forEach((value: number, xIdx: number) => {
        echartsData.push([xIdx, yIdx, value]);
        if (value < minVal) minVal = value;
        if (value > maxVal) maxVal = value;
      });
    });

    return {
      tooltip: {
        trigger: 'item',
        formatter: (params: any) => `X: ${xCenters[params.data[0]]} Å<br/>Y: ${yCenters[params.data[1]]} Å<br/>Density: ${params.data[2].toExponential(1)}`,
      },
      grid: { left: 5, right: 5, top: 5, bottom: 30, containLabel: true },
      xAxis: {
        type: 'category',
        data: xCenters,
        name: 'X (Å)',
        nameLocation: 'middle',
        nameGap: 20,
        nameTextStyle: { fontSize: 11, fontWeight: 500 },
        axisLabel: { interval: Math.floor(xCenters.length / 6), fontSize: 9 },
      },
      yAxis: {
        type: 'category',
        data: yCenters,
        name: 'Y (Å)',
        nameLocation: 'middle',
        nameGap: 25,
        nameTextStyle: { fontSize: 11, fontWeight: 500 },
        axisLabel: { interval: Math.floor(yCenters.length / 6), fontSize: 9 },
      },
      visualMap: {
        min: minVal,
        max: maxVal,
        calculable: true,
        orient: 'horizontal',
        left: 'center',
        bottom: -15,
        itemHeight: 250,
        itemWidth: 8,
        textStyle: { fontSize: 9, textAlign: "center" },
        formatter: (v: number) => v.toExponential(1),
        handleLabel: { show: true, formatter: (v: number) => v.toExponential(1) },
        inRange: { color: colorScales.viridis },
      },
      series: [{
        type: 'heatmap',
        data: echartsData,
        progressive: 0,
        progressiveThreshold: Infinity,
        animation: false,
      }],
    };
  }, [densityMap]);

  return (
    <Box sx={{ width: '100%', aspectRatio: { xs: '1 / 1.1', sm: '1 / 1' }, minHeight: { xs: 280, sm: 350 }, position: 'relative' }}>
      {loading && (
        <Box sx={{ position: 'absolute', top: 0, left: 0, right: 0, bottom: 0, backgroundColor: 'rgba(255, 255, 255, 0.7)', display: 'flex', alignItems: 'center', justifyContent: 'center', zIndex: 10 }}>
          <CircularProgress size={40} />
        </Box>
      )}
      {!densityMap && !loading && (
        <Box sx={{ height: '100%', p: 2 }}>
          <Skeleton variant="rectangular" height="100%" sx={{ borderRadius: 1 }} />
        </Box>
      )}
      {densityMap && !option && (
        <Alert severity="info" sx={{ py: 0.5 }}>No data for frames {startFrame}-{endFrame}</Alert>
      )}
      {option && (
        <ReactECharts option={option} style={{ width: '100%', height: '100%' }} opts={{ renderer: 'canvas' }} notMerge={true} />
      )}
    </Box>
  );
});

export default DensityMapChart;
