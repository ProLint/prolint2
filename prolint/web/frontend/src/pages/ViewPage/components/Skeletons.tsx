/** Skeleton loading components for ViewPage. */

import { Box, Skeleton, Grid, Card, Paper } from '@mui/material';

export const ChartSkeleton = ({ height = 300 }: { height?: number }) => (
  <Box sx={{ height, p: 2 }}>
    <Box sx={{ display: 'flex', justifyContent: 'space-between', mb: 2 }}>
      <Skeleton variant="text" width={100} height={20} />
      <Skeleton variant="text" width={60} height={20} />
    </Box>
    <Skeleton variant="rectangular" height={height - 60} sx={{ borderRadius: 1 }} />
    <Box sx={{ display: 'flex', justifyContent: 'space-between', mt: 1 }}>
      <Skeleton variant="text" width={40} height={16} />
      <Skeleton variant="text" width={40} height={16} />
      <Skeleton variant="text" width={40} height={16} />
    </Box>
  </Box>
);

export const HeatmapSkeleton = ({ height = 300, rows = 8 }: { height?: number; rows?: number }) => (
  <Box sx={{ height, p: 1 }}>
    <Box sx={{ display: 'flex', gap: 0.5, height: '100%' }}>
      {/* Y-axis labels */}
      <Box sx={{ display: 'flex', flexDirection: 'column', justifyContent: 'space-around', width: 40 }}>
        {Array.from({ length: Math.min(rows, 6) }).map((_, i) => (
          <Skeleton key={i} variant="text" width={35} height={14} />
        ))}
      </Box>
      {/* Main heatmap area */}
      <Box sx={{ flex: 1, display: 'flex', flexDirection: 'column', gap: 0.5 }}>
        {Array.from({ length: Math.min(rows, 6) }).map((_, i) => (
          <Skeleton key={i} variant="rectangular" height={`${100 / Math.min(rows, 6)}%`} sx={{ borderRadius: 0.5 }} />
        ))}
      </Box>
    </Box>
  </Box>
);

export const LogoPlotSkeleton = () => (
  <Box sx={{ height: 200, p: 2 }}>
    <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 0.5 }}>
      {Array.from({ length: 60 }).map((_, i) => (
        <Skeleton key={i} variant="rectangular" width={12} height={20 + Math.random() * 30} sx={{ borderRadius: 0.5 }} />
      ))}
    </Box>
    <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 0.5, mt: 1 }}>
      {Array.from({ length: 40 }).map((_, i) => (
        <Skeleton key={i} variant="rectangular" width={12} height={15 + Math.random() * 25} sx={{ borderRadius: 0.5 }} />
      ))}
    </Box>
  </Box>
);

export const MetricsSkeleton = ({ count = 6 }: { count?: number }) => (
  <Grid container spacing={1}>
    {Array.from({ length: count }).map((_, i) => (
      <Grid item xs={4} sm={2} key={i}>
        <Card variant="outlined" sx={{ p: 1, textAlign: 'center' }}>
          <Skeleton variant="text" width="60%" height={12} sx={{ mx: 'auto', mb: 0.5 }} />
          <Skeleton variant="text" width="80%" height={20} sx={{ mx: 'auto' }} />
        </Card>
      </Grid>
    ))}
  </Grid>
);

export const PageSkeleton = () => (
  <Box sx={{ p: 2 }}>
    {/* Header skeleton */}
    <Box sx={{ display: 'flex', gap: 1, mb: 2 }}>
      <Skeleton variant="rounded" width={80} height={24} />
      <Skeleton variant="rounded" width={100} height={24} />
      <Skeleton variant="rounded" width={60} height={24} />
    </Box>
    {/* Sticky controls skeleton */}
    <Paper sx={{ p: 1.5, mb: 2 }}>
      <Box sx={{ display: 'flex', gap: 2, alignItems: 'center', flexWrap: 'wrap' }}>
        <Skeleton variant="rectangular" width={100} height={32} sx={{ borderRadius: 1 }} />
        <Skeleton variant="rectangular" width={80} height={32} sx={{ borderRadius: 1 }} />
        <Skeleton variant="rectangular" width={80} height={32} sx={{ borderRadius: 1 }} />
        <Skeleton variant="rectangular" width={120} height={32} sx={{ borderRadius: 1 }} />
        <Skeleton variant="rectangular" width={80} height={36} sx={{ borderRadius: 1, ml: 'auto' }} />
      </Box>
    </Paper>
    {/* Content sections skeleton */}
    <Grid container spacing={2}>
      <Grid item xs={12} md={6}>
        <Paper sx={{ p: 2 }}>
          <Skeleton variant="text" width={150} height={24} sx={{ mb: 2 }} />
          <Skeleton variant="rectangular" height={250} sx={{ borderRadius: 1 }} />
        </Paper>
      </Grid>
      <Grid item xs={12} md={6}>
        <Paper sx={{ p: 2 }}>
          <Skeleton variant="text" width={120} height={24} sx={{ mb: 2 }} />
          <Box sx={{ display: 'flex', gap: 1, flexWrap: 'wrap' }}>
            {Array.from({ length: 4 }).map((_, i) => (
              <Skeleton key={i} variant="rounded" width={80} height={28} />
            ))}
          </Box>
          <Skeleton variant="rectangular" height={180} sx={{ borderRadius: 1, mt: 2 }} />
        </Paper>
      </Grid>
    </Grid>
  </Box>
);
