import { motion, AnimatePresence } from "framer-motion";
import { useState, useEffect } from "react";

interface TransitionAnimationProps {
  isAnimating: boolean;
  onComplete: () => void;
}

export const TransitionAnimation = ({ isAnimating, onComplete }: TransitionAnimationProps) => {
  const [phase, setPhase] = useState<"explode" | "reassemble" | "complete">("explode");

  useEffect(() => {
    if (isAnimating) {
      setPhase("explode");
      const timer1 = setTimeout(() => setPhase("reassemble"), 800);
      const timer2 = setTimeout(() => {
        setPhase("complete");
        onComplete();
      }, 1600);

      return () => {
        clearTimeout(timer1);
        clearTimeout(timer2);
      };
    }
  }, [isAnimating, onComplete]);

  const particles = Array.from({ length: 20 }, (_, i) => ({
    id: i,
    angle: (i / 20) * Math.PI * 2,
    distance: 100 + Math.random() * 150,
    size: 8 + Math.random() * 16,
    delay: Math.random() * 0.2,
  }));

  return (
    <AnimatePresence>
      {isAnimating && (
        <motion.div
          initial={{ opacity: 0 }}
          animate={{ opacity: 1 }}
          exit={{ opacity: 0 }}
          className="fixed inset-0 z-50 bg-background flex items-center justify-center overflow-hidden"
        >
          {/* Background gradient pulse */}
          <motion.div
            className="absolute inset-0"
            animate={{
              background: [
                "radial-gradient(circle at center, hsl(var(--primary) / 0.1) 0%, transparent 50%)",
                "radial-gradient(circle at center, hsl(var(--primary) / 0.3) 0%, transparent 70%)",
                "radial-gradient(circle at center, hsl(var(--primary) / 0.1) 0%, transparent 50%)",
              ],
            }}
            transition={{ duration: 1.6, times: [0, 0.5, 1] }}
          />

          {/* Central molecule */}
          <div className="relative">
            {/* Particles */}
            {particles.map((particle) => (
              <motion.div
                key={particle.id}
                className="absolute rounded-full"
                style={{
                  width: particle.size,
                  height: particle.size,
                  background: particle.id % 2 === 0 
                    ? 'radial-gradient(circle at 30% 30%, hsl(var(--neon-pink)), hsl(var(--neon-purple)))' 
                    : 'radial-gradient(circle at 30% 30%, hsl(var(--neon-purple)), hsl(var(--primary)))',
                  boxShadow: '0 0 20px hsl(var(--primary) / 0.6)',
                  left: '50%',
                  top: '50%',
                }}
                initial={{ x: 0, y: 0, scale: 1 }}
                animate={
                  phase === "explode"
                    ? {
                        x: Math.cos(particle.angle) * particle.distance,
                        y: Math.sin(particle.angle) * particle.distance,
                        scale: [1, 1.5, 0.8],
                        opacity: [1, 1, 0.5],
                      }
                    : phase === "reassemble"
                    ? {
                        x: 0,
                        y: 0,
                        scale: [0.8, 1.2, 1],
                        opacity: [0.5, 1, 1],
                      }
                    : {}
                }
                transition={{
                  duration: 0.8,
                  delay: particle.delay,
                  ease: phase === "explode" ? "easeOut" : "easeIn",
                }}
              />
            ))}

            {/* Central glow */}
            <motion.div
              className="w-32 h-32 rounded-full relative"
              animate={
                phase === "explode"
                  ? { scale: [1, 0.5], opacity: [1, 0] }
                  : phase === "reassemble"
                  ? { scale: [0.5, 1.2, 1], opacity: [0, 1, 1] }
                  : {}
              }
              transition={{ duration: 0.8 }}
            >
              <div 
                className="absolute inset-0 rounded-full"
                style={{
                  background: 'radial-gradient(circle, hsl(var(--primary)), transparent)',
                  filter: 'blur(20px)',
                }}
              />
            </motion.div>
          </div>

          {/* Text overlay */}
          <motion.div
            className="absolute bottom-20 text-center"
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ delay: 0.5 }}
          >
            <h2 className="font-display text-2xl font-bold gradient-text">
              {phase === "explode" ? "Analyzing..." : phase === "reassemble" ? "Initializing..." : "Welcome"}
            </h2>
          </motion.div>
        </motion.div>
      )}
    </AnimatePresence>
  );
};
