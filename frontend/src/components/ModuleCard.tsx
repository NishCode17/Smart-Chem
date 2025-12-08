import { motion } from "framer-motion";
import { useNavigate } from "react-router-dom";
import { LucideIcon } from "lucide-react";

interface ModuleCardProps {
  title: string;
  description: string;
  icon: LucideIcon;
  route: string;
  gradient?: string;
  delay?: number;
}

export const ModuleCard = ({
  title,
  description,
  icon: Icon,
  route,
  delay = 0,
}: ModuleCardProps) => {
  const navigate = useNavigate();

  return (
    <motion.div
      initial={{ opacity: 0, y: 30 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.6, delay }}
      onClick={() => navigate(route)}
      className="module-card group cursor-pointer relative overflow-hidden"
    >
      {/* Icon container - Flat & Professional */}
      <div className="relative w-14 h-14 mb-5 rounded-lg bg-primary/10 flex items-center justify-center transition-colors group-hover:bg-primary/20">
        <Icon className="w-7 h-7 text-primary/90" />
      </div>

      {/* Content */}
      <h3 className="font-display text-lg font-semibold text-foreground mb-3 transition-colors duration-200 group-hover:text-primary">
        {title}
      </h3>
      <p className="text-muted-foreground text-sm leading-relaxed">
        {description}
      </p>
    </motion.div>
  );
};
