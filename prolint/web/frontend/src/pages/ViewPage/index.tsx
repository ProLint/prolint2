/** ViewPage - Main component for viewing interaction analysis results. */

import { useMemo, useCallback } from 'react';
import {
  Grid,
  Paper,
  Typography,
  Box,
  Alert,
  CircularProgress,
  TextField,
  Button,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  Divider,
  ToggleButtonGroup,
  ToggleButton,
  Chip,
  Collapse,
  IconButton,
} from '@mui/material';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import ExpandLessIcon from '@mui/icons-material/ExpandLess';
import TouchAppIcon from '@mui/icons-material/TouchApp';
import BubbleChartIcon from '@mui/icons-material/BubbleChart';
import UploadFileIcon from '@mui/icons-material/UploadFile';

import MolecularViewer3D from '../../components/MolecularViewer3D';
import InteractionViewer3D from '../../components/InteractionViewer3D';

import { useViewPageData } from './hooks/useViewPageData';
import { getUnitLabel } from './utils';
import { PageSkeleton } from './components/Skeletons';
import DensityMapChart from './components/DensityMapChart';
import LogoPlotSection from './components/LogoPlotSection';
import { SharedContactsMatrix, SharedContactsNetworkView } from './components/SharedContactsSection';
import { TimeSeriesHeatmap, DatabaseContactsHeatmap } from './components/TimeSeriesSection';
import {
  DistanceTimePlot,
  SpatialPositionsPlot,
  AtomDistanceHeatmap,
  KineticsMetrics,
  SurvivalCurve,
  ResidenceDistribution,
} from './components/PairAnalysisSection';

export default function ViewPage() {
  const {
    resultId,
    setResultId,
    loading,
    error,
    setError,
    interactionData,
    startFrame,
    setStartFrame,
    endFrame,
    setEndFrame,
    densityMap,
    selectedDatabaseType,
    structureError,
    threeDProjection,
    threeDProjectionLoading,
    densityMapLoading,
    selectedMetric,
    logoPlotData,
    logoPlotLoading,
    sharedContactsData,
    sharedContactsLoading,
    timeSeriesData,
    timeSeriesLoading,
    selectedTimeSeriesResidue,
    residueTimeSeriesData,
    residueTimeSeriesLoading,
    selectedDatabaseResidue,
    distanceTimeData,
    distanceTimeLoading,
    distanceType,
    setDistanceType,
    pairFrame,
    pairStructure,
    pairStructureLoading,
    atomDistancesData,
    atomDistancesLoading,
    kineticsData,
    kineticsLoading,
    kineticsMode,
    selectedResidues,
    setSelectedResidues,
    expandedSections,
    compositionData,
    networkData,
    logoPlotNormalization,
    residuesByChain,
    overviewRef,
    residueSelectionRef,
    timeSeriesRef,
    pairAnalysisRef,
    isDraggingRef,
    visualSelectionRef,
    handleLoadData,
    handleDownloadStructure,
    handleFrameRangeChange,
    handleMetricChange,
    handleResidueMouseDown,
    handleResidueMouseEnter,
    handleResidueMouseUp,
    handleResidueMouseLeave,
    handleDatabaseTypeClick,
    handleTimeSeriesRowClick,
    handleDatabaseContactsRowClick,
    handlePairFrameSelect,
    handleKineticsModeChange,
    toggleSection,
    navigateToSection,
    clearTimeSeriesSelection,
  } = useViewPageData();

  // Frame range helpers
  const normalizeBy = interactionData?.params?.normalize_by || 'counts';
  const units = interactionData?.params?.units || 'ns';
  const normFactor = interactionData?.params?.norm_factor || 1;
  const showAsTime = normalizeBy === 'actual_time';
  const unitLabel = showAsTime ? getUnitLabel(units) : '';
  const frameToTime = (frame: number) => frame * normFactor;
  const timeToFrame = (time: number) => Math.round(time / normFactor);
  const minFrame = interactionData?.frame_range?.start ?? 0;
  const maxFrame = interactionData?.frame_range?.end ?? (interactionData?.universe?.n_frames ? interactionData.universe.n_frames - 1 : 100);

  const handleFrameInputChange = useCallback((setter: (v: number) => void) => (e: React.FocusEvent<HTMLInputElement>) => {
    const val = parseFloat(e.target.value);
    if (!isNaN(val)) {
      const frame = showAsTime ? timeToFrame(val) : Math.round(val);
      setter(Math.max(minFrame, Math.min(maxFrame, frame)));
    }
  }, [showAsTime, timeToFrame, minFrame, maxFrame]);

  // Memoized composition bar
  const compositionBar = useMemo(() => {
    if (compositionData.length === 0) return <Typography color="text.secondary">No data</Typography>;

    // Show dropdown when more than 5 database types
    if (compositionData.length > 5) {
      return (
        <FormControl size="small" sx={{ minWidth: 200 }}>
          <Select
            value={selectedDatabaseType || ''}
            onChange={(e) => handleDatabaseTypeClick(e.target.value as string)}
            displayEmpty
            renderValue={(value) => {
              if (!value) return <Typography color="text.secondary">Select type...</Typography>;
              const item = compositionData.find(d => d.name === value);
              return item ? `${item.name} (${item.percentage.toFixed(0)}%)` : value;
            }}
          >
            {compositionData.map((item) => (
              <MenuItem key={item.name} value={item.name}>
                <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, width: '100%' }}>
                  <Box sx={{ width: 12, height: 12, borderRadius: '50%', backgroundColor: item.color, flexShrink: 0 }} />
                  <Typography sx={{ flex: 1 }}>{item.name}</Typography>
                  <Typography variant="caption" color="text.secondary">{item.percentage.toFixed(0)}%</Typography>
                </Box>
              </MenuItem>
            ))}
          </Select>
        </FormControl>
      );
    }

    // Show chips when 5 or fewer database types
    return (
      <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 1 }}>
        {compositionData.map((item) => (
          <Chip
            key={item.name}
            label={`${item.name} (${item.percentage.toFixed(0)}%)`}
            onClick={() => handleDatabaseTypeClick(item.name)}
            sx={{
              backgroundColor: selectedDatabaseType === item.name ? item.color : 'transparent',
              color: selectedDatabaseType === item.name ? '#fff' : 'text.primary',
              border: `2px solid ${item.color}`,
              fontWeight: selectedDatabaseType === item.name ? 600 : 400,
              opacity: selectedDatabaseType === null || selectedDatabaseType === item.name ? 1 : 0.5,
              '&:hover': { backgroundColor: item.color, color: '#fff', opacity: 1 },
            }}
          />
        ))}
      </Box>
    );
  }, [compositionData, selectedDatabaseType, handleDatabaseTypeClick]);

  // Memoized selection summary
  const selectionSummary = useMemo(() => {
    if (!selectedDatabaseType) return null;
    interface SelectionItem {
      label: string;
      value: string;
    }
    const items: SelectionItem[] = [
      { label: 'Type', value: selectedDatabaseType },
      ...(selectedResidues.length > 0 ? [{ label: 'Residues', value: `${selectedResidues.length} selected` }] : []),
      ...(selectedTimeSeriesResidue ? [{ label: 'Query', value: `Res ${selectedTimeSeriesResidue}` }] : []),
      ...(selectedDatabaseResidue ? [{ label: 'Target', value: `${selectedDatabaseType} ${selectedDatabaseResidue}` }] : []),
    ];

    return (
      <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, px: 1.5, py: 0.75, bgcolor: 'grey.100', borderRadius: 1, flexWrap: 'wrap' }}>
        <Typography variant="caption" color="text.secondary" sx={{ mr: 1 }}>Selection:</Typography>
        {items.map((item, i) => (
          <Box key={i} sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
            {i > 0 && <Typography variant="caption" color="text.secondary">→</Typography>}
            <Typography variant="caption" sx={{ fontWeight: 500 }}>{item.value}</Typography>
          </Box>
        ))}
        <Button size="small" sx={{ ml: 'auto', py: 0, minWidth: 'auto', fontSize: '0.7rem' }} onClick={() => handleDatabaseTypeClick(selectedDatabaseType)}>Clear</Button>
      </Box>
    );
  }, [selectedDatabaseType, selectedResidues.length, selectedTimeSeriesResidue, selectedDatabaseResidue, handleDatabaseTypeClick]);

  return (
    <Box sx={{ height: { xs: 'auto', md: 'calc(100vh - 100px)' }, minHeight: { xs: 'calc(100vh - 100px)', md: 'auto' }, display: 'flex', flexDirection: 'column', overflow: { xs: 'visible', md: 'hidden' }, '@supports (min-height: 100dvh)': { minHeight: { xs: 'calc(100dvh - 100px)', md: 'auto' } } }}>
      {/* Header Row */}
      <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', mb: 1.5, flexShrink: 0, flexWrap: 'wrap', gap: { xs: 0.5, sm: 1 } }}>
        <Box sx={{ display: { xs: 'none', sm: 'block' } }}>
          <Typography variant="h5" sx={{ fontWeight: 600, fontSize: { xs: '1.25rem', sm: '1.5rem' } }}>Interaction Analysis</Typography>
          {!interactionData && <Typography variant="body2" color="text.secondary">Load results to explore biomolecular interactions</Typography>}
        </Box>

        {selectedDatabaseType && (
          <Box sx={{ display: 'flex', gap: 0.5, alignItems: 'center', flexWrap: 'wrap' }}>
            <Chip label="Overview" size="small" variant={expandedSections.overview ? 'filled' : 'outlined'} onClick={() => navigateToSection('overview')} sx={{ height: { xs: 28, sm: 24 }, fontSize: { xs: '0.65rem', sm: '0.7rem' } }} />
            <Chip label="Residues" size="small" variant={expandedSections.residueSelection ? 'filled' : 'outlined'} onClick={() => navigateToSection('residueSelection')} sx={{ height: { xs: 28, sm: 24 }, fontSize: { xs: '0.65rem', sm: '0.7rem' } }} />
            <Chip label="Time Series" size="small" variant={expandedSections.timeSeries ? 'filled' : 'outlined'} onClick={() => navigateToSection('timeSeries')} disabled={selectedResidues.length === 0} sx={{ height: { xs: 28, sm: 24 }, fontSize: { xs: '0.65rem', sm: '0.7rem' } }} />
            <Chip label="Pair Analysis" size="small" variant={expandedSections.pairAnalysis ? 'filled' : 'outlined'} onClick={() => navigateToSection('pairAnalysis')} disabled={!selectedTimeSeriesResidue || !selectedDatabaseResidue} sx={{ height: { xs: 28, sm: 24 }, fontSize: { xs: '0.65rem', sm: '0.7rem' } }} />
          </Box>
        )}

        {selectionSummary}
      </Box>

      {error && <Alert severity="error" sx={{ mb: 1.5, flexShrink: 0 }} onClose={() => setError(null)}>{error}</Alert>}

      {/* Initial Load Section */}
      {!interactionData && !loading && (
        <Paper sx={{ p: 4, textAlign: 'center', flex: 1, display: 'flex', flexDirection: 'column', justifyContent: 'center', alignItems: 'center' }}>
          <UploadFileIcon sx={{ fontSize: 48, color: 'grey.400', mb: 2 }} />
          <Typography variant="h6" sx={{ mb: 1, fontWeight: 500 }}>Load Computation Results</Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 3, maxWidth: 400 }}>Enter the result ID from your computation to explore biomolecular interactions</Typography>
          <Box sx={{ display: 'flex', gap: 1, alignItems: 'center', width: '100%', maxWidth: 400 }}>
            <TextField fullWidth size="small" label="Result ID" value={resultId} onChange={(e) => setResultId(e.target.value)} placeholder="e.g., abc123-def456" />
            <Button variant="contained" onClick={() => handleLoadData()} disabled={!resultId} sx={{ minWidth: 80 }}>Load</Button>
          </Box>
        </Paper>
      )}

      {loading && <PageSkeleton />}

      {interactionData && !loading && (
        <Box sx={{ flex: 1, overflow: 'auto', pr: 1 }}>
          {/* Controls Bar */}
          <Paper sx={{ p: 1.5, mb: 1.5, position: 'sticky', top: 0, zIndex: 10, backgroundColor: 'background.paper', borderBottom: '1px solid', borderColor: 'divider' }} elevation={1}>
            <Grid container spacing={1.5} alignItems="center">
              <Grid item xs={12} md={5}>
                <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                  <Typography variant="caption" color="text.secondary" sx={{ whiteSpace: 'nowrap' }}>Select type:</Typography>
                  {compositionBar}
                </Box>
              </Grid>

              {selectedDatabaseType && (
                <>
                  <Grid item xs={6} sm={3} md={1.5}>
                    <TextField label={showAsTime ? `Start (${unitLabel})` : "Start"} type="number" size="small" fullWidth defaultValue={showAsTime ? frameToTime(startFrame) : startFrame} key={`start-${showAsTime}-${startFrame}`} onBlur={handleFrameInputChange(setStartFrame)} inputProps={{ min: showAsTime ? frameToTime(minFrame) : minFrame, max: showAsTime ? frameToTime(maxFrame) : maxFrame, step: 'any' }} />
                  </Grid>
                  <Grid item xs={6} sm={3} md={1.5}>
                    <TextField label={showAsTime ? `End (${unitLabel})` : "End"} type="number" size="small" fullWidth defaultValue={showAsTime ? frameToTime(endFrame) : endFrame} key={`end-${showAsTime}-${endFrame}`} onBlur={handleFrameInputChange(setEndFrame)} inputProps={{ min: showAsTime ? frameToTime(minFrame) : minFrame, max: showAsTime ? frameToTime(maxFrame) : maxFrame, step: 'any' }} />
                  </Grid>
                  <Grid item xs={6} sm={3} md={2}>
                    <FormControl fullWidth size="small">
                      <InputLabel>Metric</InputLabel>
                      <Select value={selectedMetric} label="Metric" onChange={(e) => handleMetricChange(e.target.value)}>
                        <MenuItem value="occupancy">Occupancy</MenuItem>
                        <MenuItem value="mean">Mean</MenuItem>
                        <MenuItem value="max">Max</MenuItem>
                      </Select>
                    </FormControl>
                  </Grid>
                  <Grid item xs={6} sm={3} md={1.5}>
                    <Button variant="contained" fullWidth size="small" onClick={handleFrameRangeChange} disabled={densityMapLoading || threeDProjectionLoading}>
                      {densityMapLoading || threeDProjectionLoading ? <CircularProgress size={18} color="inherit" /> : 'Update'}
                    </Button>
                  </Grid>
                </>
              )}
            </Grid>
          </Paper>

          {/* Empty State */}
          {!selectedDatabaseType && (
            <Paper sx={{ p: 4, textAlign: 'center', mt: 2 }}>
              <TouchAppIcon sx={{ fontSize: 48, color: 'primary.main', mb: 2, opacity: 0.7 }} />
              <Typography variant="h6" sx={{ mb: 1, fontWeight: 500 }}>Select a Database Residue Type</Typography>
              <Typography variant="body2" color="text.secondary" sx={{ maxWidth: 400, mx: 'auto' }}>Click on a database residue type above to explore density maps, 3D structures, time series, kinetics, and more.</Typography>
            </Paper>
          )}

          {/* Visualizations */}
          {selectedDatabaseType && (
            <>
              {/* Section 1: Overview */}
              <Paper ref={overviewRef} sx={{ mb: 1.5, overflow: 'hidden' }}>
                <Box onClick={() => toggleSection('overview')} sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', p: 1, cursor: 'pointer', bgcolor: 'grey.50', borderBottom: expandedSections.overview ? '1px solid' : 'none', borderColor: 'divider', '&:hover': { bgcolor: 'grey.100' } }}>
                  <Typography variant="subtitle2" sx={{ fontWeight: 600 }}>System Overview</Typography>
                  <IconButton size="small" sx={{ p: 0.25 }}>{expandedSections.overview ? <ExpandLessIcon fontSize="small" /> : <ExpandMoreIcon fontSize="small" />}</IconButton>
                </Box>
                <Collapse in={expandedSections.overview}>
                  <Box sx={{ p: { xs: 0.5, sm: 1 } }}>
                    <Grid container spacing={{ xs: 1, sm: 1.5 }}>
                      <Grid item xs={12} sm={6}>
                        <Typography variant="caption" sx={{ fontWeight: 600, display: 'block', mb: 0, textAlign: 'center' }}>2D Density Map</Typography>
                        <DensityMapChart densityMap={densityMap} loading={densityMapLoading} startFrame={startFrame} endFrame={endFrame} />
                      </Grid>
                      <Grid item xs={12} sm={6}>
                        <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', mb: 0.5, position: 'relative' }}>
                          <Typography variant="caption" sx={{ fontWeight: 600 }}>3D Structure</Typography>
                          <Button size="small" variant="text" onClick={(e) => { e.stopPropagation(); handleDownloadStructure(); }} sx={{ position: 'absolute', right: 0, fontSize: '0.65rem' }}>Download</Button>
                        </Box>
                        <MolecularViewer3D pdbData={threeDProjection} metric={selectedMetric} loading={threeDProjectionLoading} error={structureError} />
                      </Grid>
                    </Grid>
                  </Box>
                </Collapse>
              </Paper>

              {/* Section 2: Residue Selection */}
              <Paper ref={residueSelectionRef} sx={{ mb: 1.5, overflow: 'hidden' }}>
                <Box onClick={() => toggleSection('residueSelection')} sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', p: 1, cursor: 'pointer', bgcolor: 'grey.50', borderBottom: expandedSections.residueSelection ? '1px solid' : 'none', borderColor: 'divider', '&:hover': { bgcolor: 'grey.100' } }}>
                  <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                    <Typography variant="subtitle2" sx={{ fontWeight: 600 }}>Residue Selection</Typography>
                    {selectedResidues.length > 0 && <Chip label={`${selectedResidues.length} selected`} size="small" color="primary" sx={{ height: 20, fontSize: '0.7rem' }} />}
                  </Box>
                  <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                    {selectedResidues.length > 0 && <Button size="small" variant="text" onClick={(e) => { e.stopPropagation(); setSelectedResidues([]); }} sx={{ fontSize: '0.7rem', py: 0 }}>Clear</Button>}
                    <IconButton size="small" sx={{ p: 0.25 }}>{expandedSections.residueSelection ? <ExpandLessIcon fontSize="small" /> : <ExpandMoreIcon fontSize="small" />}</IconButton>
                  </Box>
                </Box>
                <Collapse in={expandedSections.residueSelection}>
                  <Box sx={{ p: 1 }}>
                    <Typography variant="caption" color="text.secondary" sx={{ display: 'block', mb: 0.5 }}>Click/drag to select residues</Typography>
                    <LogoPlotSection
                      logoPlotData={logoPlotData}
                      loading={logoPlotLoading}
                      selectedMetric={selectedMetric}
                      selectedResidues={selectedResidues}
                      logoPlotNormalization={logoPlotNormalization}
                      residuesByChain={residuesByChain}
                      isDraggingRef={isDraggingRef}
                      visualSelectionRef={visualSelectionRef}
                      onMouseDown={handleResidueMouseDown}
                      onMouseEnter={handleResidueMouseEnter}
                      onMouseUp={handleResidueMouseUp}
                      onMouseLeave={handleResidueMouseLeave}
                    />
                    <Grid container spacing={{ xs: 1, sm: 1.5 }} sx={{ mt: 1 }}>
                      <Grid item xs={12} sm={6}>
                        <Typography variant="caption" sx={{ fontWeight: 600, display: 'block', mb: 0, textAlign: 'center' }}>Shared Contacts Matrix</Typography>
                        <Box sx={{ width: '100%', aspectRatio: { xs: '1 / 1.2', sm: '1 / 1' }, minHeight: { xs: 220, sm: 280 } }}>
                          <SharedContactsMatrix data={sharedContactsData} loading={sharedContactsLoading} />
                        </Box>
                      </Grid>
                      <Grid item xs={12} sm={6}>
                        <Typography variant="caption" sx={{ fontWeight: 600, display: 'block', mb: 0.5, textAlign: 'center' }}>Shared Contacts Network</Typography>
                        <Box sx={{ minHeight: { xs: 220, sm: 280 }, display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
                          {selectedResidues.length === 0 ? (
                            <Box sx={{ textAlign: 'center', p: 2 }}>
                              <BubbleChartIcon sx={{ fontSize: 36, color: 'grey.400', mb: 1 }} />
                              <Typography variant="body2" color="text.secondary">Select residues from the logo plot above</Typography>
                            </Box>
                          ) : (
                            <SharedContactsNetworkView data={sharedContactsData} loading={sharedContactsLoading} networkData={networkData} selectedResiduesCount={selectedResidues.length} />
                          )}
                        </Box>
                      </Grid>
                    </Grid>
                  </Box>
                </Collapse>
              </Paper>

              {/* Section 3: Time Series */}
              {selectedResidues.length > 0 && (
                <Paper ref={timeSeriesRef} sx={{ mb: 1.5, overflow: 'hidden' }}>
                  <Box onClick={() => toggleSection('timeSeries')} sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', p: 1, cursor: 'pointer', bgcolor: 'grey.50', borderBottom: expandedSections.timeSeries ? '1px solid' : 'none', borderColor: 'divider', '&:hover': { bgcolor: 'grey.100' } }}>
                    <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                      <Typography variant="subtitle2" sx={{ fontWeight: 600 }}>Time Series</Typography>
                      {selectedTimeSeriesResidue && <Chip label={`Res ${selectedTimeSeriesResidue}`} size="small" color="secondary" sx={{ height: 20, fontSize: '0.7rem' }} />}
                    </Box>
                    <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                      {selectedTimeSeriesResidue && <Button size="small" variant="text" onClick={(e) => { e.stopPropagation(); clearTimeSeriesSelection(); }} sx={{ fontSize: '0.7rem', py: 0 }}>Clear</Button>}
                      <IconButton size="small" sx={{ p: 0.25 }}>{expandedSections.timeSeries ? <ExpandLessIcon fontSize="small" /> : <ExpandMoreIcon fontSize="small" />}</IconButton>
                    </Box>
                  </Box>
                  <Collapse in={expandedSections.timeSeries}>
                    <Box sx={{ p: 1 }}>
                      <Typography variant="caption" color="text.secondary" sx={{ display: 'block', mb: 1 }}>Click row to explore</Typography>
                      <TimeSeriesHeatmap data={timeSeriesData} loading={timeSeriesLoading} interactionData={interactionData} onRowClick={handleTimeSeriesRowClick} />
                      {selectedTimeSeriesResidue && (
                        <>
                          <Divider sx={{ my: 1.5 }} />
                          <Typography variant="caption" sx={{ fontWeight: 600, display: 'block', mb: 0.5 }}>Database Contacts for Residue {selectedTimeSeriesResidue}</Typography>
                          <Typography variant="caption" color="text.secondary" sx={{ display: 'block', mb: 1 }}>Click row for pair analysis</Typography>
                          <DatabaseContactsHeatmap data={residueTimeSeriesData} loading={residueTimeSeriesLoading} interactionData={interactionData} selectedDatabaseType={selectedDatabaseType} onRowClick={handleDatabaseContactsRowClick} />
                        </>
                      )}
                    </Box>
                  </Collapse>
                </Paper>
              )}

              {/* Section 4: Pair Analysis */}
              {selectedTimeSeriesResidue && selectedDatabaseResidue && (
                <Paper ref={pairAnalysisRef} sx={{ mb: 1.5, overflow: 'hidden' }}>
                  <Box onClick={() => toggleSection('pairAnalysis')} sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', p: 1, cursor: 'pointer', bgcolor: 'grey.50', borderBottom: expandedSections.pairAnalysis ? '1px solid' : 'none', borderColor: 'divider', '&:hover': { bgcolor: 'grey.100' } }}>
                    <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                      <Typography variant="subtitle2" sx={{ fontWeight: 600 }}>Pair Analysis</Typography>
                      <Chip label={`Res ${selectedTimeSeriesResidue} ↔ ${selectedDatabaseType} ${selectedDatabaseResidue}`} size="small" color="success" sx={{ height: 20, fontSize: '0.65rem' }} />
                    </Box>
                    <IconButton size="small" sx={{ p: 0.25 }}>{expandedSections.pairAnalysis ? <ExpandLessIcon fontSize="small" /> : <ExpandMoreIcon fontSize="small" />}</IconButton>
                  </Box>
                  <Collapse in={expandedSections.pairAnalysis}>
                    <Box sx={{ p: 1 }}>
                      {distanceTimeData && (
                        <>
                          <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 0.5 }}>
                            <Typography variant="caption" sx={{ fontWeight: 600 }}>Distance vs Time</Typography>
                            <ToggleButtonGroup value={distanceType} exclusive onChange={(_, v) => v && setDistanceType(v)} size="small">
                              <ToggleButton value="min" sx={{ px: 1, py: 0, fontSize: '0.65rem' }}>Min</ToggleButton>
                              <ToggleButton value="com" sx={{ px: 1, py: 0, fontSize: '0.65rem' }}>COM</ToggleButton>
                            </ToggleButtonGroup>
                          </Box>
                          <DistanceTimePlot data={distanceTimeData} loading={distanceTimeLoading} distanceType={distanceType} pairFrame={pairFrame} interactionData={interactionData} onFrameSelect={handlePairFrameSelect} />

                          <Divider sx={{ my: 1.5 }} />
                          <Grid container spacing={{ xs: 1, sm: 1.5 }}>
                            <Grid item xs={12} sm={6}>
                              <Typography variant="caption" sx={{ fontWeight: 600, display: 'block', mb: 0.5, textAlign: 'center' }}>2D Positions (Frame {distanceTimeData.frames[pairFrame]})</Typography>
                              <SpatialPositionsPlot distanceTimeData={distanceTimeData} densityMap={densityMap} pairFrame={pairFrame} selectedTimeSeriesResidue={selectedTimeSeriesResidue} selectedDatabaseResidue={selectedDatabaseResidue} selectedDatabaseType={selectedDatabaseType} />
                            </Grid>
                            <Grid item xs={12} sm={6}>
                              <Typography variant="caption" sx={{ fontWeight: 600, display: 'block', mb: 0.5, textAlign: 'center' }}>3D View (Frame {distanceTimeData.frames[pairFrame]})</Typography>
                              <InteractionViewer3D pdbData={pairStructure} queryResidue={selectedTimeSeriesResidue!} databaseResidue={selectedDatabaseResidue} loading={pairStructureLoading} />
                            </Grid>
                          </Grid>

                          <Divider sx={{ my: 1.5 }} />
                          <Grid container spacing={{ xs: 1, sm: 1.5 }} sx={{ alignItems: 'stretch' }}>
                            <Grid item xs={12} sm={6} sx={{ display: 'flex', flexDirection: 'column' }}>
                              <Typography variant="caption" sx={{ fontWeight: 600, display: 'block', mb: 0, textAlign: 'center' }}>Atom-Atom Distances</Typography>
                              <Box sx={{ flex: 1, minHeight: { xs: 200, sm: 250 } }}>
                                <AtomDistanceHeatmap data={atomDistancesData} loading={atomDistancesLoading} selectedTimeSeriesResidue={selectedTimeSeriesResidue} selectedDatabaseResidue={selectedDatabaseResidue} />
                              </Box>
                            </Grid>
                            <Grid item xs={12} sm={6} sx={{ display: 'flex', flexDirection: 'column' }}>
                              <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'center', gap: 1, mb: 0.5 }}>
                                <Typography variant="caption" sx={{ fontWeight: 600, mb: 0.5 }}>Kinetics</Typography>
                                <ToggleButtonGroup value={kineticsMode} exclusive onChange={(_, value) => value && handleKineticsModeChange(value)} size="small">
                                  <ToggleButton value="individual" sx={{ px: 0.75, py: 0, fontSize: '0.6rem' }}>Single</ToggleButton>
                                  <ToggleButton value="accumulated" sx={{ px: 0.75, py: 0, fontSize: '0.6rem' }}>All {selectedDatabaseType}</ToggleButton>
                                </ToggleButtonGroup>
                              </Box>
                              <KineticsMetrics data={kineticsData} loading={kineticsLoading} interactionData={interactionData} />
                              <Divider sx={{ my: 1 }} />
                              <Grid container spacing={1} sx={{ flex: 1, minHeight: { xs: 160, sm: 200 } }}>
                                <Grid item xs={6} sx={{ display: 'flex', flexDirection: 'column' }}>
                                  <Typography variant="caption" sx={{ fontWeight: 500, display: 'block', textAlign: 'center', fontSize: '0.65rem' }}>Survival</Typography>
                                  <Box sx={{ flex: 1, minHeight: { xs: 120, sm: 150 } }}>
                                    <SurvivalCurve data={kineticsData} interactionData={interactionData} />
                                  </Box>
                                </Grid>
                                <Grid item xs={6} sx={{ display: 'flex', flexDirection: 'column' }}>
                                  <Typography variant="caption" sx={{ fontWeight: 500, display: 'block', textAlign: 'center', fontSize: '0.65rem' }}>Distribution</Typography>
                                  <Box sx={{ flex: 1, minHeight: { xs: 120, sm: 150 } }}>
                                    <ResidenceDistribution data={kineticsData} interactionData={interactionData} />
                                  </Box>
                                </Grid>
                              </Grid>
                            </Grid>
                          </Grid>
                        </>
                      )}
                    </Box>
                  </Collapse>
                </Paper>
              )}
            </>
          )}
        </Box>
      )}
    </Box>
  );
}
