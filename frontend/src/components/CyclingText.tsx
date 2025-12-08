import { motion, AnimatePresence } from "framer-motion";
import { useEffect, useState } from "react";

const words = ["Research", "Innovate", "Design", "Discover"];

export const CyclingText = () => {
  const [currentIndex, setCurrentIndex] = useState(0);

  useEffect(() => {
    const interval = setInterval(() => {
      setCurrentIndex((prev) => (prev + 1) % words.length);
    }, 2000);
    return () => clearInterval(interval);
  }, []);

  return (
    <div className="h-20 flex items-center justify-start overflow-hidden">
      <AnimatePresence mode="wait">
        <motion.span
          key={currentIndex}
          className="text-5xl md:text-6xl lg:text-7xl font-display font-bold gradient-text"
          initial={{ y: 50, opacity: 0 }}
          animate={{ y: 0, opacity: 1 }}
          exit={{ y: -50, opacity: 0 }}
          transition={{ duration: 0.5, ease: "easeOut" }}
        >
          {words[currentIndex]}
        </motion.span>
      </AnimatePresence>
    </div>
  );
};
