/** Shared contacts network visualization using Cytoscape.js. */

import { useEffect, useRef, useState, useCallback } from 'react';
import { Box, Typography, Alert, Chip, Stack, Button } from '@mui/material';
import ArrowBackIcon from '@mui/icons-material/ArrowBack';
import cytoscape from 'cytoscape';
import { visualizationTheme } from '../theme/visualizationTheme';
import { interpolateColorScale } from '../utils/colorUtils';

interface Node {
  id: string;
  label: string;
  type: string;
  residue_id: number;
  restype: string;
}

interface Edge {
  source: string;
  target: string;
  value: number;
  shared_db_residues?: number[];
}

interface SharedContactsNetworkProps {
  nodes: Node[];
  edges: Edge[];
}

export default function SharedContactsNetwork({ nodes, edges }: SharedContactsNetworkProps) {
  const cyRef = useRef<HTMLDivElement>(null);
  const cyInstanceRef = useRef<cytoscape.Core | null>(null);
  const [focusedNode, setFocusedNode] = useState<string | null>(null);
  const [focusedNodeData, setFocusedNodeData] = useState<Node | null>(null);
  const [networkStats, setNetworkStats] = useState<{
    totalEdges: number;
    maxShared: number;
    avgShared: number;
  } | null>(null);

  // Back to full network handler
  const handleBackToFull = useCallback(() => {
    const cy = cyInstanceRef.current;
    if (!cy) return;

    setFocusedNode(null);
    setFocusedNodeData(null);

    // Show all elements
    cy.elements().removeClass('hidden focused neighbor selected');
    cy.elements().style('display', 'element');

    // Re-layout and fit
    cy.layout({
      name: 'circle',
      animate: true,
      animationDuration: 300,
      padding: 20,
      fit: true,
    }).run();
  }, []);

  // Zoom to node handler
  const handleZoomToNode = useCallback((nodeId: string) => {
    const cy = cyInstanceRef.current;
    if (!cy) return;

    const node = cy.getElementById(nodeId);
    if (!node || node.length === 0) return;

    // Get connected elements
    const connectedEdges = node.connectedEdges();
    const neighbors = node.neighborhood('node');
    const visibleElements = node.union(connectedEdges).union(neighbors);

    // Hide non-connected elements
    cy.elements().style('display', 'none');
    visibleElements.style('display', 'element');

    // Style the focused node and neighbors
    cy.elements().removeClass('focused neighbor selected');
    node.addClass('focused');
    neighbors.addClass('neighbor');
    connectedEdges.addClass('selected');

    // Layout only visible elements in a concentric pattern (focused in center)
    visibleElements.layout({
      name: 'concentric',
      animate: true,
      animationDuration: 300,
      fit: true,
      padding: 20,
      concentric: (n: cytoscape.NodeSingular) => n.id() === nodeId ? 2 : 1,
      levelWidth: () => 1,
      minNodeSpacing: 30,
    }).run();

    setFocusedNode(nodeId);
    setFocusedNodeData(nodes.find(n => n.id === nodeId) || null);
  }, [nodes]);

  // Initialize network
  useEffect(() => {
    if (!cyRef.current || !nodes || nodes.length === 0) return;

    // Calculate statistics
    const edgeValues = edges.map(e => e.value);
    const stats = {
      totalEdges: edges.length,
      maxShared: Math.max(...edgeValues, 1),
      avgShared: edgeValues.reduce((a, b) => a + b, 0) / Math.max(edges.length, 1),
    };
    setNetworkStats(stats);

    // Convert data to Cytoscape format
    const elements = [
      ...nodes.map(node => ({
        data: {
          id: node.id,
          label: node.label,
          residue_id: node.residue_id,
          restype: node.restype,
        },
      })),
      ...edges.map((edge, idx) => ({
        data: {
          id: `edge-${idx}`,
          source: edge.source,
          target: edge.target,
          value: edge.value,
          width: 1 + (edge.value / stats.maxShared) * 4,
          color: interpolateColorScale(edge.value, 1, stats.maxShared, visualizationTheme.colorscales.sharedContacts),
        },
      })),
    ];

    // Create Cytoscape instance
    const cy = cytoscape({
      container: cyRef.current,
      elements,
      style: visualizationTheme.cytoscapeStyle as cytoscape.StylesheetStyle[],
      layout: {
        name: 'circle',
        animate: true,
        animationDuration: 500,
        padding: 30,
        fit: true,
        avoidOverlap: true,
      },
      wheelSensitivity: 0.2,
      minZoom: 0.1,
      maxZoom: 3,
    });

    // Fit and center the graph after layout
    cy.on('layoutstop', () => {
      cy.resize();
      cy.fit(undefined, 10);
      cy.center();
    });

    // Initial fit after a short delay to ensure container is rendered
    setTimeout(() => {
      cy.resize();
      cy.fit(undefined, 10);
      cy.center();
    }, 200);

    // Handle window resize
    const handleResize = () => {
      if (cy) {
        cy.resize();
        cy.fit(undefined, 10);
        cy.center();
      }
    };
    window.addEventListener('resize', handleResize);

    cyInstanceRef.current = cy;

    return () => {
      window.removeEventListener('resize', handleResize);
      if (cyInstanceRef.current) {
        cyInstanceRef.current.destroy();
      }
    };
  }, [nodes, edges]);

  // Separate effect for tap handlers to avoid stale closures
  useEffect(() => {
    const cy = cyInstanceRef.current;
    if (!cy) return;

    const handleNodeTap = (event: cytoscape.EventObject) => {
      const node = event.target;
      const nodeId = node.data('id');

      // If already focused on this node, go back to full
      if (focusedNode === nodeId) {
        handleBackToFull();
      } else {
        handleZoomToNode(nodeId);
      }
    };

    const handleBackgroundTap = (event: cytoscape.EventObject) => {
      if (event.target === cy && focusedNode) {
        handleBackToFull();
      }
    };

    cy.on('tap', 'node', handleNodeTap);
    cy.on('tap', handleBackgroundTap);

    return () => {
      cy.off('tap', 'node', handleNodeTap);
      cy.off('tap', handleBackgroundTap);
    };
  }, [focusedNode, handleBackToFull, handleZoomToNode]);

  // Get connection count for focused node
  const focusedNodeConnections = focusedNode
    ? edges.filter(e => e.source === focusedNode || e.target === focusedNode).length
    : 0;

  if (!nodes || nodes.length === 0) {
    return (
      <Alert severity="info">
        No shared contacts found. Select more residues above.
      </Alert>
    );
  }

  // Generate accessible description
  const networkDescription = networkStats
    ? `Network graph showing ${nodes.length} residues with ${networkStats.totalEdges} connections. ${focusedNode ? `Focused on ${focusedNodeData?.label || focusedNode}.` : 'Click a node to focus on it.'}`
    : `Network graph with ${nodes.length} residues`;

  return (
    <Box sx={{ overflow: 'hidden' }}>
      {/* Single Cytoscape container for both views */}
      <Box
        role="img"
        aria-label={networkDescription}
        sx={{
          position: 'relative',
          border: `1px solid ${visualizationTheme.colors.border}`,
          borderRadius: 2,
          overflow: 'hidden',
          backgroundColor: visualizationTheme.colors.surface,
        }}
      >
        <div
          ref={cyRef}
          style={{
            width: '100%',
            height: `${Math.min(visualizationTheme.plotDimensions.networkHeight, 450)}px`,
          }}
          aria-hidden="true"
        />
      </Box>

      <Typography variant="caption" color="text.secondary" sx={{ mt: 1, display: 'block', textAlign: 'center' }}>
        {focusedNode ? (
          <>Click node or background to return</>
        ) : (
          <>Click node to focus. Scroll to zoom, drag to pan.</>
        )}
      </Typography>

      {/* Stats/info at the bottom */}
      <Box sx={{ mt: 1.5, display: 'flex', justifyContent: 'center' }}>
        {focusedNode ? (
          <Stack direction="row" spacing={1} flexWrap="wrap" useFlexGap sx={{ justifyContent: 'center', alignItems: 'center' }}>
            <Button
              size="small"
              startIcon={<ArrowBackIcon />}
              onClick={handleBackToFull}
              variant="outlined"
            >
              Back to full network
            </Button>
            <Typography variant="body2" color="text.secondary">
              <strong>{focusedNodeData?.label}</strong> — {focusedNodeConnections} connection{focusedNodeConnections !== 1 ? 's' : ''}
            </Typography>
          </Stack>
        ) : (
          networkStats && (
            <Stack direction="row" spacing={1} flexWrap="wrap" useFlexGap sx={{ justifyContent: 'center' }}>
              <Chip
                label={`${nodes.length} Residue${nodes.length !== 1 ? 's' : ''}`}
                color="primary"
                variant="outlined"
                size="small"
              />
              <Chip
                label={`${networkStats.totalEdges} Connection${networkStats.totalEdges !== 1 ? 's' : ''}`}
                color="primary"
                variant="outlined"
                size="small"
              />
              {networkStats.totalEdges > 0 && (
                <>
                  <Chip
                    label={`Max: ${networkStats.maxShared} shared`}
                    color="secondary"
                    variant="outlined"
                    size="small"
                  />
                  <Chip
                    label={`Avg: ${networkStats.avgShared.toFixed(1)} shared`}
                    color="secondary"
                    variant="outlined"
                    size="small"
                  />
                  <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
                    <Box
                      sx={{
                        width: 32,
                        height: 3,
                        background: `linear-gradient(to right,
                          ${visualizationTheme.colorscales.sharedContacts[0][1]},
                          ${visualizationTheme.colorscales.sharedContacts[1][1]},
                          ${visualizationTheme.colorscales.sharedContacts[2][1]})`,
                        borderRadius: 1,
                        flexShrink: 0,
                      }}
                    />
                    <Typography variant="caption">Low → High Shared</Typography>
                  </Box>
                </>
              )}
            </Stack>
          )
        )}
      </Box>
    </Box>
  );
}

