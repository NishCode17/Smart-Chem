import { motion } from "framer-motion";
import { useNavigate } from "react-router-dom";
import { api } from "@/lib/api";
import { ModuleCard } from "@/components/ModuleCard";
import { MoleculeBackground } from "@/components/MoleculeBackground";
import {
  FlaskConical,
  Microscope,
  MessageSquare,
  Library,
  LogOut,
  Settings,
  User
} from "lucide-react";
import { Button } from "@/components/ui/button";
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuSeparator,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu";

const Dashboard = () => {
  const navigate = useNavigate();

  const modules = [
    {
      title: "Design Studio",
      description: "Generate and optimize novel molecules with AI-powered tools. Create drug candidates from scratch or enhance existing leads.",
      icon: FlaskConical,
      route: "/design-studio",
    },
    {
      title: "Virtual Lab",
      description: "Analyze molecules in 3D, view ADMET properties, and get AI-driven insights from Dr. SmartChem.",
      icon: Microscope,
      route: "/virtual-lab",
    },
    {
      title: "Smart Assistant",
      description: "Chat with Dr. SmartChem for chemistry insights, literature search, and drug design guidance.",
      icon: MessageSquare,
      route: "/assistant",
    },
    {
      title: "Library",
      description: "Organize your molecules in projects, search your collection, and export data for collaboration.",
      icon: Library,
      route: "/library",
    },
  ];

  return (
    <div className="min-h-screen bg-background relative overflow-hidden">
      <MoleculeBackground />

      {/* Header */}
      <header className="relative z-10 border-b border-border/50 backdrop-blur-xl bg-background/50">
        <div className="container mx-auto px-6 py-4 flex items-center justify-between">
          <motion.div
            initial={{ opacity: 0, x: -20 }}
            animate={{ opacity: 1, x: 0 }}
            className="flex items-center gap-3"
          >
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-primary to-secondary flex items-center justify-center">
              <FlaskConical className="w-5 h-5 text-primary-foreground" />
            </div>
            <span className="font-display text-xl font-bold text-foreground">
              SmartChem <span className="gradient-text">AI</span>
            </span>
          </motion.div>

          <motion.div
            initial={{ opacity: 0, x: 20 }}
            animate={{ opacity: 1, x: 0 }}
            className="flex items-center gap-4"
          >
            <Button variant="ghost" size="icon" className="rounded-full">
              <Settings className="w-5 h-5" />
            </Button>

            <DropdownMenu>
              <DropdownMenuTrigger asChild>
                <Button variant="glass" size="icon" className="rounded-full">
                  <User className="w-5 h-5" />
                </Button>
              </DropdownMenuTrigger>
              <DropdownMenuContent align="end" className="w-48">
                <DropdownMenuItem>
                  <User className="w-4 h-4 mr-2" />
                  Profile
                </DropdownMenuItem>
                <DropdownMenuItem>
                  <Settings className="w-4 h-4 mr-2" />
                  Settings
                </DropdownMenuItem>
                <DropdownMenuSeparator />
                <DropdownMenuItem onClick={() => api.logout()}>
                  <LogOut className="w-4 h-4 mr-2" />
                  Sign Out
                </DropdownMenuItem>
              </DropdownMenuContent>
            </DropdownMenu>
          </motion.div>
        </div>
      </header>

      {/* Main Content */}
      <main className="relative z-10 container mx-auto px-6 py-12">
        {/* Welcome Section */}
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          className="text-center mb-16"
        >
          <h1 className="font-display text-4xl md:text-5xl font-bold text-foreground mb-4">
            Welcome to <span className="gradient-text">SmartChem AI</span>
          </h1>
          <p className="text-lg text-muted-foreground max-w-2xl mx-auto">
            Your AI-powered drug discovery platform. Generate molecules, analyze properties,
            and accelerate your research with cutting-edge machine learning.
          </p>
        </motion.div>

        {/* Module Grid */}
        <div className="grid grid-cols-1 md:grid-cols-2 gap-6 max-w-4xl mx-auto">
          {modules.map((module, index) => (
            <ModuleCard
              key={module.title}
              {...module}
              delay={0.1 + index * 0.1}
            />
          ))}
        </div>

        {/* Quick Stats */}
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ delay: 0.6 }}
          className="mt-16 grid grid-cols-2 md:grid-cols-4 gap-4 max-w-4xl mx-auto"
        >
          {[
            { label: "Molecules Generated", value: "1,247" },
            { label: "Lead Candidates", value: "23" },
            { label: "Projects", value: "8" },
            { label: "AI Analyses", value: "456" },
          ].map((stat, i) => (
            <div key={i} className="glass-card p-4 text-center">
              <p className="font-display text-2xl font-bold gradient-text">{stat.value}</p>
              <p className="text-xs text-muted-foreground mt-1">{stat.label}</p>
            </div>
          ))}
        </motion.div>
      </main>
    </div>
  );
};

export default Dashboard;
