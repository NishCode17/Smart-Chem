import { motion } from "framer-motion";
import { PageHeader } from "@/components/PageHeader";
import { ChatInterface } from "@/components/ChatInterface";
import { Card } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { 
  FlaskConical, 
  FileText, 
  Lightbulb, 
  Search,
  Sparkles
} from "lucide-react";

const Assistant = () => {
  const quickActions = [
    {
      icon: FlaskConical,
      label: "Analyze a molecule",
      prompt: "Can you analyze this molecule and tell me about its drug-likeness?",
    },
    {
      icon: Search,
      label: "Search literature",
      prompt: "Find recent papers about EGFR inhibitors in cancer treatment",
    },
    {
      icon: Lightbulb,
      label: "Optimization ideas",
      prompt: "Suggest ways to improve the bioavailability of my lead compound",
    },
    {
      icon: FileText,
      label: "Explain a concept",
      prompt: "Explain the Lipinski Rule of 5 and why it matters",
    },
  ];

  return (
    <div className="min-h-screen bg-background p-6">
      <div className="container mx-auto max-w-5xl">
        <PageHeader
          title="Smart Assistant"
          description="Chat with Dr. SmartChem for AI-powered chemistry insights"
        />

        <div className="grid grid-cols-1 lg:grid-cols-4 gap-6">
          {/* Quick Actions Sidebar */}
          <motion.div
            initial={{ opacity: 0, x: -20 }}
            animate={{ opacity: 1, x: 0 }}
            className="lg:col-span-1 space-y-3"
          >
            <h3 className="font-display text-sm font-semibold text-foreground mb-4 flex items-center gap-2">
              <Sparkles className="w-4 h-4 text-primary" />
              Quick Actions
            </h3>
            
            {quickActions.map((action, i) => (
              <motion.div
                key={i}
                initial={{ opacity: 0, x: -10 }}
                animate={{ opacity: 1, x: 0 }}
                transition={{ delay: i * 0.1 }}
              >
                <Button
                  variant="outline"
                  className="w-full justify-start text-left h-auto py-3 px-4"
                >
                  <action.icon className="w-4 h-4 mr-3 text-primary flex-shrink-0" />
                  <span className="text-sm">{action.label}</span>
                </Button>
              </motion.div>
            ))}

            {/* Tips Card */}
            <Card className="p-4 mt-6">
              <h4 className="font-display text-sm font-semibold text-foreground mb-2">
                Pro Tips
              </h4>
              <ul className="space-y-2 text-xs text-muted-foreground">
                <li className="flex items-start gap-2">
                  <span className="text-primary">•</span>
                  Paste SMILES strings directly for analysis
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-primary">•</span>
                  Ask for ADMET predictions
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-primary">•</span>
                  Request synthesis routes
                </li>
                <li className="flex items-start gap-2">
                  <span className="text-primary">•</span>
                  Compare multiple molecules
                </li>
              </ul>
            </Card>
          </motion.div>

          {/* Chat Interface */}
          <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            className="lg:col-span-3"
          >
            <Card className="h-[calc(100vh-220px)] min-h-[500px] overflow-hidden">
              <ChatInterface />
            </Card>
          </motion.div>
        </div>
      </div>
    </div>
  );
};

export default Assistant;
