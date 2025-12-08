import { motion } from "framer-motion";
import { useEffect, useState } from "react";

interface Particle {
  id: number;
  x: number;
  y: number;
  size: number;
  duration: number;
  delay: number;
}

interface MoleculeNode {
  id: number;
  x: number;
  y: number;
  size: number;
  color: string;
}

interface MoleculeBond {
  from: number;
  to: number;
}

export const MoleculeBackground = () => {
  const [particles, setParticles] = useState<Particle[]>([]);
  const [molecules, setMolecules] = useState<{ nodes: MoleculeNode[]; bonds: MoleculeBond[] }[]>([]);

  useEffect(() => {
    // Generate particles
    const newParticles: Particle[] = Array.from({ length: 50 }, (_, i) => ({
      id: i,
      x: Math.random() * 100,
      y: Math.random() * 100,
      size: Math.random() * 4 + 2,
      duration: Math.random() * 10 + 10,
      delay: Math.random() * 5,
    }));
    setParticles(newParticles);

    // Generate molecule structures
    const moleculeConfigs = [
      { centerX: 20, centerY: 30, scale: 1 },
      { centerX: 70, centerY: 60, scale: 0.8 },
      { centerX: 40, centerY: 80, scale: 0.6 },
      { centerX: 85, centerY: 20, scale: 0.7 },
    ];

    const newMolecules = moleculeConfigs.map((config, molIndex) => {
      const nodeCount = Math.floor(Math.random() * 3) + 4;
      const nodes: MoleculeNode[] = [];
      const bonds: MoleculeBond[] = [];

      for (let i = 0; i < nodeCount; i++) {
        const angle = (i / nodeCount) * Math.PI * 2;
        const radius = 5 * config.scale;
        nodes.push({
          id: i,
          x: config.centerX + Math.cos(angle) * radius,
          y: config.centerY + Math.sin(angle) * radius,
          size: (Math.random() * 1 + 0.5) * config.scale,
          color: i % 2 === 0 ? "var(--neon-purple)" : "var(--neon-pink)",
        });

        if (i > 0) {
          bonds.push({ from: i - 1, to: i });
        }
      }
      bonds.push({ from: nodeCount - 1, to: 0 });

      return { nodes, bonds };
    });

    setMolecules(newMolecules);
  }, []);

  return (
    <div className="absolute inset-0 overflow-hidden pointer-events-none">
      {/* Gradient overlay */}
      <div className="absolute inset-0 bg-gradient-to-br from-primary/5 via-transparent to-secondary/5" />
      
      {/* Floating particles */}
      {particles.map((particle) => (
        <motion.div
          key={particle.id}
          className="absolute rounded-full bg-primary/30"
          style={{
            left: `${particle.x}%`,
            top: `${particle.y}%`,
            width: particle.size,
            height: particle.size,
          }}
          animate={{
            y: [0, -30, 0],
            x: [0, 10, -10, 0],
            opacity: [0.3, 0.7, 0.3],
          }}
          transition={{
            duration: particle.duration,
            delay: particle.delay,
            repeat: Infinity,
            ease: "easeInOut",
          }}
        />
      ))}

      {/* Molecule structures */}
      <svg className="absolute inset-0 w-full h-full">
        {molecules.map((molecule, molIndex) => (
          <motion.g
            key={molIndex}
            animate={{
              rotate: [0, 360],
            }}
            transition={{
              duration: 30 + molIndex * 10,
              repeat: Infinity,
              ease: "linear",
            }}
            style={{
              transformOrigin: `${molecule.nodes[0]?.x}% ${molecule.nodes[0]?.y}%`,
            }}
          >
            {/* Bonds */}
            {molecule.bonds.map((bond, bondIndex) => {
              const from = molecule.nodes[bond.from];
              const to = molecule.nodes[bond.to];
              return (
                <motion.line
                  key={bondIndex}
                  x1={`${from.x}%`}
                  y1={`${from.y}%`}
                  x2={`${to.x}%`}
                  y2={`${to.y}%`}
                  stroke="hsl(var(--primary) / 0.3)"
                  strokeWidth="1"
                  initial={{ pathLength: 0 }}
                  animate={{ pathLength: 1 }}
                  transition={{ duration: 2, delay: molIndex * 0.5 }}
                />
              );
            })}
            {/* Nodes */}
            {molecule.nodes.map((node) => (
              <motion.circle
                key={node.id}
                cx={`${node.x}%`}
                cy={`${node.y}%`}
                r={node.size * 5}
                fill={`hsl(${node.color} / 0.6)`}
                animate={{
                  scale: [1, 1.2, 1],
                  opacity: [0.4, 0.8, 0.4],
                }}
                transition={{
                  duration: 3,
                  delay: node.id * 0.2,
                  repeat: Infinity,
                  ease: "easeInOut",
                }}
              />
            ))}
          </motion.g>
        ))}
      </svg>
    </div>
  );
};
