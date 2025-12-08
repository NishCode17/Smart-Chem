import { useEffect, useRef, useState } from "react";
import { motion } from "framer-motion";
import { Maximize2, RotateCcw, ZoomIn, ZoomOut } from "lucide-react";
import { Button } from "@/components/ui/button";
// @ts-ignore
import * as $3Dmol from "3dmol";
import { api } from "@/lib/api";

interface MoleculeViewerProps {
  smiles?: string;
  pdbData?: string;
  className?: string;
}

export const MoleculeViewer = ({ smiles, pdbData, className }: MoleculeViewerProps) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const viewerRef = useRef<any>(null);

  useEffect(() => {
    if (!containerRef.current) return;

    // Initialize Viewer
    const element = containerRef.current;
    const config = { backgroundColor: "#1a1a2e" }; // Match --deep-space color roughly
    const viewer = $3Dmol.createViewer(element, config);
    viewerRef.current = viewer;

    const loadMolecule = async () => {
      try {
        setIsLoading(true);
        setError(null);
        viewer.clear();

        let molData = pdbData;

        if (smiles && !molData) {
          // Fetch 3D structure from backend
          const response = await api.get3DStructure(smiles);
          molData = response.mol_block;
        }

        if (molData) {
          viewer.addModel(molData, "mol");
          viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.25 } });
          viewer.zoomTo();
          viewer.render();
        } else {
          // If nothing to load, we can leave it empty or show placeholder
          // (Logic handled by conditional rendering below if needed, but here we just clear)
        }
      } catch (err) {
        console.error("3D Viewer Error:", err);
        setError("Failed to load 3D structure");
      } finally {
        setIsLoading(false);
      }
    };

    if (smiles || pdbData) {
      loadMolecule();
    }

    return () => {
      // Cleanup if needed
    };
  }, [smiles, pdbData]);

  const handleZoomIn = () => {
    if (viewerRef.current) {
      viewerRef.current.zoom(1.2, 1000);
    }
  };

  const handleZoomOut = () => {
    if (viewerRef.current) {
      viewerRef.current.zoom(0.8, 1000);
    }
  };

  const handleReset = () => {
    if (viewerRef.current) {
      viewerRef.current.zoomTo();
    }
  };

  // Placeholder 3D visualization (Only shown if NO smiles provided initially)
  const renderPlaceholderMolecule = () => (
    <div className="absolute inset-0 flex items-center justify-center">
      <motion.div
        className="relative"
        animate={{ rotateY: 360 }}
        transition={{ duration: 20, repeat: Infinity, ease: "linear" }}
      >
        {/* Central atom */}
        <div className="w-12 h-12 rounded-full molecule-node absolute left-1/2 top-1/2 -translate-x-1/2 -translate-y-1/2" />

        {/* Surrounding atoms */}
        {[0, 72, 144, 216, 288].map((angle, i) => (
          <motion.div
            key={i}
            className="absolute"
            style={{
              left: `calc(50% + ${Math.cos((angle * Math.PI) / 180) * 60}px)`,
              top: `calc(50% + ${Math.sin((angle * Math.PI) / 180) * 60}px)`,
            }}
            animate={{ scale: [1, 1.2, 1] }}
            transition={{ duration: 2, delay: i * 0.2, repeat: Infinity }}
          >
            <div
              className="w-8 h-8 rounded-full -translate-x-1/2 -translate-y-1/2"
              style={{
                background: i % 2 === 0
                  ? 'radial-gradient(circle at 30% 30%, hsl(var(--neon-pink)), hsl(var(--neon-purple)))'
                  : 'radial-gradient(circle at 30% 30%, hsl(var(--neon-blue)), hsl(var(--neon-purple)))',
                boxShadow: '0 0 15px hsl(var(--neon-purple) / 0.5)'
              }}
            />
            {/* Bond line */}
            <svg className="absolute inset-0 w-full h-full pointer-events-none" style={{ overflow: 'visible' }}>
              <line
                x1="0"
                y1="0"
                x2={-Math.cos((angle * Math.PI) / 180) * 45}
                y2={-Math.sin((angle * Math.PI) / 180) * 45}
                stroke="hsl(var(--primary) / 0.5)"
                strokeWidth="2"
              />
            </svg>
          </motion.div>
        ))}
      </motion.div>
    </div>
  );

  return (
    <div className={`relative bg-background rounded-xl overflow-hidden ${className}`}>
      {/* Toolbar */}
      <div className="absolute top-3 right-3 z-10 flex gap-2">
        <Button variant="glass" size="icon" className="w-8 h-8" onClick={handleZoomIn}>
          <ZoomIn className="w-4 h-4" />
        </Button>
        <Button variant="glass" size="icon" className="w-8 h-8" onClick={handleZoomOut}>
          <ZoomOut className="w-4 h-4" />
        </Button>
        <Button variant="glass" size="icon" className="w-8 h-8" onClick={handleReset}>
          <RotateCcw className="w-4 h-4" />
        </Button>
        <Button variant="glass" size="icon" className="w-8 h-8">
          <Maximize2 className="w-4 h-4" />
        </Button>
      </div>

      {/* Viewer container */}
      <div
        ref={containerRef}
        className="w-full h-full min-h-[300px] relative"
      >
        {isLoading && (
          <div className="absolute inset-0 flex items-center justify-center z-20 bg-background/50 backdrop-blur-sm">
            <motion.div
              animate={{ rotate: 360 }}
              transition={{ duration: 2, repeat: Infinity, ease: "linear" }}
              className="w-12 h-12 border-2 border-primary border-t-transparent rounded-full"
            />
          </div>
        )}

        {error ? (
          <div className="absolute inset-0 flex items-center justify-center text-destructive">
            {error}
          </div>
        ) : (
          (!smiles && !pdbData) && renderPlaceholderMolecule()
        )}
      </div>

      {/* SMILES display */}
      {smiles && (
        <div className="absolute bottom-3 left-3 right-3 bg-muted/80 backdrop-blur-sm rounded-lg px-3 py-2 z-10">
          <span className="font-mono text-xs text-muted-foreground truncate block">
            {smiles}
          </span>
        </div>
      )}

      {/* Grid overlay */}
      <div
        className="absolute inset-0 pointer-events-none opacity-10"
        style={{
          backgroundImage: 'radial-gradient(circle at 1px 1px, hsl(var(--primary)) 1px, transparent 0)',
          backgroundSize: '40px 40px'
        }}
      />
    </div>
  );
};
