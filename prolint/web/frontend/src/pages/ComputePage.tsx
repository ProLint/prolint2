/** ComputePage - Upload files and run interaction computations. */

import React, { useState, useRef, useCallback, memo, useEffect } from 'react';
import {
  Typography,
  Paper,
  TextField,
  Button,
  Box,
  Alert,
  CircularProgress,
  IconButton,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  Chip,
  Divider,
  LinearProgress,
} from '@mui/material';
import CloudUploadIcon from '@mui/icons-material/CloudUpload';
import InsertDriveFileIcon from '@mui/icons-material/InsertDriveFile';
import DeleteIcon from '@mui/icons-material/Delete';
import PlayArrowIcon from '@mui/icons-material/PlayArrow';
import VisibilityIcon from '@mui/icons-material/Visibility';
import CheckCircleIcon from '@mui/icons-material/CheckCircle';
import axios from 'axios';
import { useNavigate } from 'react-router-dom';
import type { DatasetInfo, ComputeResponse, ReplicaAnalysisResult } from '../types';

// Accepted file extensions
const TOPOLOGY_EXTENSIONS = ['.gro', '.pdb', '.tpr', '.crd', '.psf', '.mol2', '.xyz'];
const TRAJECTORY_EXTENSIONS = ['.xtc', '.trr', '.dcd', '.nc', '.netcdf', '.lammpstrj'];

// Maximum number of frames allowed for analysis to avoid timeout issues
const MAX_FRAMES_LIMIT = 5000;

interface CompactFileUploadProps {
  label: string;
  file: File | null;
  onFileSelect: (file: File | null) => void;
  acceptedExtensions: string[];
  required?: boolean;
}

const CompactFileUpload = memo(function CompactFileUpload({ label, file, onFileSelect, acceptedExtensions, required }: CompactFileUploadProps) {
  const inputRef = useRef<HTMLInputElement>(null);
  const [isDragOver, setIsDragOver] = useState(false);

  const handleDragOver = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    setIsDragOver(true);
  }, []);

  const handleDragLeave = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    setIsDragOver(false);
  }, []);

  const handleDrop = useCallback((e: React.DragEvent) => {
    e.preventDefault();
    setIsDragOver(false);
    const droppedFile = e.dataTransfer.files[0];
    if (droppedFile) onFileSelect(droppedFile);
  }, [onFileSelect]);

  const formatFileSize = (bytes: number) => {
    if (bytes < 1024) return `${bytes} B`;
    if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
    return `${(bytes / (1024 * 1024)).toFixed(1)} MB`;
  };

  return (
    <Box>
      <Typography variant="caption" color="text.secondary" sx={{ mb: 0.5, display: 'block' }}>
        {label} {required && <span style={{ color: '#d32f2f' }}>*</span>}
      </Typography>
      <input
        ref={inputRef}
        type="file"
        accept={acceptedExtensions.join(',')}
        onChange={(e) => e.target.files?.[0] && onFileSelect(e.target.files[0])}
        style={{ display: 'none' }}
      />
      <Box
        onClick={() => inputRef.current?.click()}
        onDragOver={handleDragOver}
        onDragLeave={handleDragLeave}
        onDrop={handleDrop}
        sx={{
          border: '1px dashed',
          borderColor: isDragOver ? 'primary.main' : file ? 'success.main' : 'grey.300',
          borderRadius: 1,
          p: 1.5,
          cursor: 'pointer',
          backgroundColor: isDragOver ? 'action.hover' : file ? 'success.50' : 'grey.50',
          transition: 'all 0.15s ease',
          '&:hover': { borderColor: 'primary.main', backgroundColor: 'action.hover' },
          minHeight: 48,
          display: 'flex',
          alignItems: 'center',
        }}
      >
        {file ? (
          <Box sx={{ display: 'flex', alignItems: 'center', width: '100%', gap: 1 }}>
            <InsertDriveFileIcon sx={{ fontSize: 20, color: 'success.main' }} />
            <Box sx={{ flex: 1, minWidth: 0 }}>
              <Typography variant="body2" noWrap sx={{ fontWeight: 500 }}>
                {file.name}
              </Typography>
              <Typography variant="caption" color="text.secondary">
                {formatFileSize(file.size)}
              </Typography>
            </Box>
            <IconButton
              size="small"
              onClick={(e) => { e.stopPropagation(); onFileSelect(null); }}
              color="error"
            >
              <DeleteIcon fontSize="small" />
            </IconButton>
          </Box>
        ) : (
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, color: 'text.secondary' }}>
            <CloudUploadIcon sx={{ fontSize: 20 }} />
            <Typography variant="body2">
              Drop file or click to browse
            </Typography>
          </Box>
        )}
      </Box>
    </Box>
  );
});

export default function ComputePage() {
  const navigate = useNavigate();
  const [loading, setLoading] = useState(false);
  const [computing, setComputing] = useState(false);
  const [analyzing, setAnalyzing] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [dataset, setDataset] = useState<DatasetInfo | null>(null);
  const [result, setResult] = useState<ComputeResponse | null>(null);
  const [replicaAnalysis, setReplicaAnalysis] = useState<ReplicaAnalysisResult | null>(null);

  // File uploads
  const [topologyFile, setTopologyFile] = useState<File | null>(null);
  const [trajectoryFile, setTrajectoryFile] = useState<File | null>(null);

  // Replica selection (only used when repeated residue IDs are detected)
  const [selectedReplica, setSelectedReplica] = useState<string | null>(null);

  // Computation parameters
  const [formData, setFormData] = useState({
    query_selection: 'protein',
    database_selection: 'not protein',
    cutoff: 7.0,
    start: 0,
    stop: 100,
    step: 1,
    units: 'ns',
    normalize_by: 'counts',
  });

  // Clear replica analysis when query_selection changes (requires re-analysis)
  useEffect(() => {
    if (dataset) {
      setReplicaAnalysis(null);
      setSelectedReplica(null);
    }
  }, [formData.query_selection]);

  const timeUnits = [
    { value: 'fs', label: 'fs' },
    { value: 'ps', label: 'ps' },
    { value: 'ns', label: 'ns' },
    { value: 'us', label: 'μs' },
    { value: 'ms', label: 'ms' },
  ];

  // Helper for number inputs - allows empty field while typing
  const handleNumberChange = useCallback((field: string, isFloat = false) => (e: React.ChangeEvent<HTMLInputElement>) => {
    const val = e.target.value;
    if (val === '') {
      setFormData(prev => ({ ...prev, [field]: '' as unknown }));
    } else {
      const num = isFloat ? parseFloat(val) : parseInt(val, 10);
      if (!isNaN(num)) {
        setFormData(prev => ({ ...prev, [field]: num }));
      }
    }
  }, []);

  // Calculate number of frames that will be analyzed
  const calculateFramesToAnalyze = useCallback((start: number, stop: number, step: number): number => {
    if (step <= 0 || stop <= start) return 0;
    return Math.ceil((stop - start) / step);
  }, []);

  // Current number of frames to analyze
  const framesToAnalyze = calculateFramesToAnalyze(
    Number(formData.start) || 0,
    Number(formData.stop) || 0,
    Number(formData.step) || 1
  );

  // Check if frame count exceeds limit
  const exceedsFrameLimit = framesToAnalyze > MAX_FRAMES_LIMIT;

  const handleLoadFiles = async () => {
    if (!topologyFile) {
      setError('Please upload a topology file');
      return;
    }

    setLoading(true);
    setError(null);

    try {
      const formDataPayload = new FormData();
      formDataPayload.append('topology_file', topologyFile);
      if (trajectoryFile) formDataPayload.append('trajectory_file', trajectoryFile);

      const response = await axios.post('/api/datasets/upload', formDataPayload, {
        headers: { 'Content-Type': 'multipart/form-data' },
      });

      setDataset(response.data);
      setReplicaAnalysis(null); // Reset replica analysis when new files are loaded
      if (response.data.n_frames) {
        // Calculate step to limit frames to MAX_FRAMES_LIMIT
        const nFrames = response.data.n_frames;
        const defaultStep = nFrames > MAX_FRAMES_LIMIT ? Math.ceil(nFrames / MAX_FRAMES_LIMIT) : 1;
        setFormData(prev => ({ ...prev, stop: nFrames, step: defaultStep }));
      }
    } catch (err: any) {
      setError(err.response?.data?.detail || err.message || 'Failed to upload files');
    } finally {
      setLoading(false);
    }
  };

  // Helper to format error messages for user display
  const formatErrorMessage = (err: any, fallback: string): string => {
    const detail = err.response?.data?.detail || err.message || fallback;
    // Make common MDAnalysis errors more user-friendly
    if (detail.includes('Selection') && detail.includes('matched no atoms')) {
      return `Invalid selection: no atoms matched. Check your selection syntax.`;
    }
    if (detail.includes('SelectionError') || detail.includes('selection')) {
      return `Invalid selection syntax. Use MDAnalysis selection language (e.g., "protein", "resname LIP*", "resid 1-10").`;
    }
    return detail;
  };

  const handleAnalyzeReplicas = async () => {
    if (!dataset) {
      setError('Please load files first');
      return;
    }

    // Validate query selection is not empty
    if (!formData.query_selection.trim()) {
      setError('Query selection cannot be empty');
      return;
    }

    setAnalyzing(true);
    setError(null);

    try {
      const response = await axios.post('/api/dashboard/analyze-replicas', {
        dataset_id: dataset.id,
        query_selection: formData.query_selection,
      });

      setReplicaAnalysis(response.data);

      // Set default replica selection when repeated residue IDs detected
      if (response.data.n_replicas > 1 && response.data.has_repeated_resids && response.data.replica_info?.length > 0) {
        // Default to first replica when selection is required
        setSelectedReplica(response.data.replica_info[0].replica_id);
      } else {
        setSelectedReplica(null);
      }
    } catch (err: any) {
      setError(formatErrorMessage(err, 'Replica analysis failed'));
    } finally {
      setAnalyzing(false);
    }
  };

  const handleCompute = async () => {
    if (!dataset) {
      setError('Please load files first');
      return;
    }

    // Validate selections are not empty
    if (!formData.query_selection.trim()) {
      setError('Query selection cannot be empty');
      return;
    }
    if (!formData.database_selection.trim()) {
      setError('Database selection cannot be empty');
      return;
    }

    // Validate frame count does not exceed limit
    if (exceedsFrameLimit) {
      setError(`Frame count (${framesToAnalyze.toLocaleString()}) exceeds the maximum limit of ${MAX_FRAMES_LIMIT.toLocaleString()} frames. Please increase the step value or reduce the frame range to continue.`);
      return;
    }

    // Validate replica selection when repeated residue IDs detected
    const needsReplicaSelection = replicaAnalysis?.n_replicas
      && replicaAnalysis.n_replicas > 1
      && replicaAnalysis?.has_repeated_resids;

    if (needsReplicaSelection && !selectedReplica) {
      setError('Please select a replica for analysis');
      return;
    }

    setComputing(true);
    setError(null);

    try {
      const payload = {
        dataset_id: dataset.id,
        ...formData,
        selected_replica: selectedReplica,
        replica_info: replicaAnalysis?.replica_info || null,
      };
      const response = await axios.post('/api/dashboard/compute', payload);
      setResult(response.data);
    } catch (err: any) {
      setError(formatErrorMessage(err, 'Computation failed'));
    } finally {
      setComputing(false);
    }
  };

  const handleReset = () => {
    setDataset(null);
    setResult(null);
    setTopologyFile(null);
    setTrajectoryFile(null);
    setSelectedReplica(null);
    setReplicaAnalysis(null);
    setFormData({
      query_selection: 'protein',
      database_selection: 'not protein',
      cutoff: 7.0,
      start: 0,
      stop: 100,
      step: 1,
      units: 'ns',
      normalize_by: 'counts',
    });
  };

  // Check if repeated residue IDs were detected (requires replica selection)
  const hasRepeatedResidueIds = !!(replicaAnalysis?.n_replicas
    && replicaAnalysis.n_replicas > 1
    && replicaAnalysis?.has_repeated_resids === true);

  // Show replica selection when repeated residue IDs are detected
  const showReplicaSelector = hasRepeatedResidueIds
    && replicaAnalysis?.replica_info
    && replicaAnalysis.replica_info.length > 0;

  // Success state
  if (result) {
    return (
      <Box sx={{
        minHeight: { xs: 'calc(100vh - 100px)', md: 'calc(100vh - 120px)' },
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
        p: { xs: 2, sm: 3 },
        '@supports (min-height: 100dvh)': {
          minHeight: { xs: 'calc(100dvh - 100px)', md: 'calc(100dvh - 120px)' },
        },
      }}>
        <Paper sx={{ p: { xs: 3, sm: 4 }, textAlign: 'center', maxWidth: 400, width: '100%' }}>
          <CheckCircleIcon sx={{ fontSize: { xs: 48, sm: 64 }, color: 'success.main', mb: 2 }} />
          <Typography variant="h5" sx={{ fontSize: { xs: '1.25rem', sm: '1.5rem' } }} gutterBottom>
            Computation Complete
          </Typography>
          <Typography color="text.secondary" sx={{ mb: 1 }}>
            Computed in {result.computation_time.toFixed(2)}s
          </Typography>
          <Typography variant="body2" color="text.secondary" sx={{ mb: 3, wordBreak: 'break-all' }}>
            Result ID: {result.result_id}
          </Typography>
          <Box sx={{ display: 'flex', gap: 2, justifyContent: 'center', flexDirection: { xs: 'column', sm: 'row' } }}>
            <Button
              variant="contained"
              startIcon={<VisibilityIcon />}
              onClick={() => navigate(`/view?result_id=${result.result_id}`)}
              fullWidth
              sx={{ maxWidth: { sm: 'auto' } }}
            >
              View Results
            </Button>
            <Button variant="outlined" onClick={handleReset} fullWidth sx={{ maxWidth: { sm: 'auto' } }}>
              New Analysis
            </Button>
          </Box>
        </Paper>
      </Box>
    );
  }

  return (
    <Box sx={{
      minHeight: { xs: 'auto', md: 'calc(100vh - 120px)' },
      display: 'flex',
      alignItems: { xs: 'flex-start', md: 'center' },
      justifyContent: 'center',
      p: { xs: 2, sm: 3 },
      py: { xs: 3, md: 3 },
      '@supports (min-height: 100dvh)': {
        minHeight: { xs: 'auto', md: 'calc(100dvh - 120px)' },
      },
    }}>
      <Paper sx={{
        p: { xs: 2.5, sm: 4 },
        width: '100%',
      }}>
        {/* Header */}
        <Typography variant="h5" sx={{ fontWeight: 600, mb: 0.5, textAlign: 'center', fontSize: { xs: '1.25rem', sm: '1.5rem' } }}>
          Compute Interactions
        </Typography>
        <Typography variant="body2" color="text.secondary" sx={{ mb: { xs: 2, sm: 3 }, textAlign: 'center' }}>
          Upload simulation files and configure analysis
        </Typography>

        {error && (
          <Alert severity="error" sx={{ mb: { xs: 2, sm: 3 } }} onClose={() => setError(null)}>
            {error}
          </Alert>
        )}

        {/* File Uploads */}
        <Box sx={{ display: 'flex', flexDirection: { xs: 'column', sm: 'row' }, gap: 2, mb: 2 }}>
          <Box sx={{ flex: 1 }}>
            <CompactFileUpload
              label="Topology / Coordinates"
              file={topologyFile}
              onFileSelect={setTopologyFile}
              acceptedExtensions={TOPOLOGY_EXTENSIONS}
              required
            />
          </Box>
          <Box sx={{ flex: 1 }}>
            <CompactFileUpload
              label="Trajectory (optional)"
              file={trajectoryFile}
              onFileSelect={setTrajectoryFile}
              acceptedExtensions={TRAJECTORY_EXTENSIONS}
            />
          </Box>
        </Box>

        {/* Load Button & Dataset Info */}
        <Box sx={{ display: 'flex', flexDirection: { xs: 'column', sm: 'row' }, alignItems: { xs: 'stretch', sm: 'center' }, gap: 2, mb: 2 }}>
          <Button
            variant={dataset ? 'outlined' : 'contained'}
            size="small"
            startIcon={loading ? <CircularProgress size={16} color="inherit" /> : <CloudUploadIcon />}
            onClick={handleLoadFiles}
            disabled={loading || !topologyFile}
            sx={{ minWidth: 120 }}
          >
            {loading ? 'Loading...' : dataset ? 'Reload' : 'Load Files'}
          </Button>
          {dataset && (
            <Box sx={{ display: 'flex', alignItems: 'center', gap: 1, flex: 1, flexWrap: 'wrap' }}>
              <CheckCircleIcon sx={{ fontSize: 18, color: 'success.main' }} />
              <Chip label={`${dataset.n_frames} frames`} size="small" variant="outlined" />
              <Chip label={`${dataset.n_atoms?.toLocaleString()} atoms`} size="small" variant="outlined" />
              <Chip label={`${dataset.n_residues?.toLocaleString()} residues`} size="small" variant="outlined" />
            </Box>
          )}
        </Box>

        <Divider sx={{ mb: { xs: 2, sm: 3 } }} />

        {/* Selections Row */}
        <Box sx={{ display: 'flex', flexDirection: { xs: 'column', sm: 'row' }, gap: { xs: 1.5, sm: 2 }, mb: 2 }}>
          <TextField
            size="small"
            fullWidth
            label="Query Selection"
            value={formData.query_selection}
            onChange={(e) => setFormData(prev => ({ ...prev, query_selection: e.target.value }))}
            placeholder="e.g., protein"
          />
          <TextField
            size="small"
            fullWidth
            label="Database Selection"
            value={formData.database_selection}
            onChange={(e) => setFormData(prev => ({ ...prev, database_selection: e.target.value }))}
            placeholder="e.g., not protein"
          />
        </Box>

        {/* Frame Range Row */}
        <Box sx={{ display: 'flex', flexDirection: { xs: 'column', sm: 'row' }, gap: { xs: 1.5, sm: 2 }, mb: 2 }}>
          <Box sx={{ display: 'flex', gap: { xs: 1.5, sm: 2 }, flex: { sm: 2 } }}>
            <TextField
              size="small"
              fullWidth
              type="number"
              label="Start Frame"
              value={formData.start}
              onChange={handleNumberChange('start')}
              inputProps={{ min: 0 }}
            />
            <TextField
              size="small"
              fullWidth
              type="number"
              label="Stop Frame"
              value={formData.stop}
              onChange={handleNumberChange('stop')}
              inputProps={{ max: dataset?.n_frames }}
            />
          </Box>
          <TextField
            size="small"
            fullWidth
            type="number"
            label="Step"
            value={formData.step}
            onChange={handleNumberChange('step')}
            inputProps={{ min: 1 }}
            sx={{ flex: { sm: 1 } }}
          />
        </Box>

        {/* Frame count info and warning */}
        <Box sx={{ mb: 2, display: 'flex', alignItems: 'center', gap: 1, flexWrap: 'wrap' }}>
          <Typography variant="body2" color="text.secondary">
            Frames to analyze: <strong>{framesToAnalyze.toLocaleString()}</strong>
          </Typography>
          {exceedsFrameLimit && (
            <Alert severity="warning" sx={{ flex: 1, py: 0.5 }}>
              Exceeds limit of {MAX_FRAMES_LIMIT.toLocaleString()} frames. Increase step or reduce range.
            </Alert>
          )}
        </Box>

        {/* Settings Row */}
        <Box sx={{ display: 'flex', flexDirection: { xs: 'column', sm: 'row' }, gap: { xs: 1.5, sm: 2 }, mb: { xs: 2, sm: 3 } }}>
          <TextField
            size="small"
            fullWidth
            type="number"
            label="Cutoff (Å)"
            value={formData.cutoff}
            onChange={handleNumberChange('cutoff', true)}
            inputProps={{ step: 0.1, min: 0.1 }}
          />
          <FormControl size="small" fullWidth>
            <InputLabel>Normalize By</InputLabel>
            <Select
              value={formData.normalize_by}
              label="Normalize By"
              onChange={(e) => setFormData(prev => ({ ...prev, normalize_by: e.target.value }))}
            >
              <MenuItem value="counts">Frame Counts</MenuItem>
              <MenuItem value="actual_time">Actual Time</MenuItem>
            </Select>
          </FormControl>
          <FormControl size="small" fullWidth disabled={formData.normalize_by !== 'actual_time'}>
            <InputLabel>Time Units</InputLabel>
            <Select
              value={formData.units}
              label="Time Units"
              onChange={(e) => setFormData(prev => ({ ...prev, units: e.target.value }))}
            >
              {timeUnits.map((unit) => (
                <MenuItem key={unit.value} value={unit.value}>{unit.label}</MenuItem>
              ))}
            </Select>
          </FormControl>
        </Box>

        {/* Replica Analysis Section */}
        <Box sx={{ mb: 2 }}>
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 2, mb: 1.5 }}>
            <Button
              variant="outlined"
              size="small"
              onClick={handleAnalyzeReplicas}
              disabled={analyzing || !dataset}
              startIcon={analyzing ? <CircularProgress size={14} /> : null}
            >
              {analyzing ? 'Analyzing...' : 'Analyze Replicas'}
            </Button>
            {replicaAnalysis && (
              <Typography variant="body2" color="text.secondary">
                {replicaAnalysis.message}
              </Typography>
            )}
          </Box>

          {/* Replica Selection - shown when repeated residue IDs detected */}
          {showReplicaSelector && (
            <Box sx={{ p: 2, bgcolor: 'grey.50', borderRadius: 1, border: '1px solid', borderColor: 'grey.200' }}>
              <Typography variant="subtitle2" sx={{ mb: 1.5, fontWeight: 600 }}>
                Select Replica for Analysis
              </Typography>
              <Typography variant="body2" color="text.secondary" sx={{ mb: 1.5 }}>
                {replicaAnalysis?.n_replicas} replicas with repeated residue IDs detected. Select which replica to analyze.
              </Typography>
              <FormControl size="small" fullWidth>
                <InputLabel>Select Replica</InputLabel>
                <Select
                  value={selectedReplica || ''}
                  label="Select Replica"
                  onChange={(e) => setSelectedReplica(e.target.value)}
                >
                  {replicaAnalysis?.replica_info?.map((replica) => (
                    <MenuItem key={replica.replica_id} value={replica.replica_id}>
                      Replica {replica.replica_id} (residues {replica.first_resid}-{replica.last_resid}, {replica.n_residues} residues)
                    </MenuItem>
                  ))}
                </Select>
              </FormControl>
            </Box>
          )}
        </Box>

        {/* Progress Bar */}
        {computing && (
          <Box sx={{ mb: 2 }}>
            <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', mb: 1 }}>
              <Typography variant="body2" color="text.secondary">
                Computing interactions...
              </Typography>
              <CircularProgress size={16} />
            </Box>
            <LinearProgress
              sx={{
                height: 8,
                borderRadius: 1,
                backgroundColor: 'grey.200',
                '& .MuiLinearProgress-bar': {
                  borderRadius: 1,
                }
              }}
            />
            <Typography variant="caption" color="text.secondary" sx={{ mt: 0.5, display: 'block' }}>
              Processing trajectory frames. This may take a moment for large files.
            </Typography>
          </Box>
        )}

        {/* Compute Button */}
        <Button
          variant="contained"
          size="large"
          fullWidth
          startIcon={computing ? <CircularProgress size={20} color="inherit" /> : <PlayArrowIcon />}
          onClick={handleCompute}
          disabled={computing || !dataset || !replicaAnalysis || (hasRepeatedResidueIds && !selectedReplica) || exceedsFrameLimit}
          sx={{ py: { xs: 1.25, sm: 1.5 } }}
        >
          {computing ? 'Computing...' : 'Compute Interactions'}
        </Button>
        {!dataset && (
          <Typography variant="caption" color="text.secondary" sx={{ mt: 1.5, display: 'block', textAlign: 'center' }}>
            Load files to continue
          </Typography>
        )}
        {dataset && !replicaAnalysis && (
          <Typography variant="caption" color="text.secondary" sx={{ mt: 1.5, display: 'block', textAlign: 'center' }}>
            Analyze replicas to continue
          </Typography>
        )}
        {replicaAnalysis && hasRepeatedResidueIds && !selectedReplica && (
          <Typography variant="caption" color="text.secondary" sx={{ mt: 1.5, display: 'block', textAlign: 'center' }}>
            Select a replica to enable computation
          </Typography>
        )}
        {replicaAnalysis && exceedsFrameLimit && (
          <Typography variant="caption" color="error" sx={{ mt: 1.5, display: 'block', textAlign: 'center' }}>
            Reduce frames to analyze (max {MAX_FRAMES_LIMIT.toLocaleString()})
          </Typography>
        )}
      </Paper>
    </Box>
  );
}
