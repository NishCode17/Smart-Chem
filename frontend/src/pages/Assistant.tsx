import { motion } from "framer-motion";
import { PageHeader } from "@/components/PageHeader";
import { ChatInterface } from "@/components/ChatInterface";
import { Card } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import {
  FlaskConical,
  BrainCircuit,
  Lightbulb,
  ShieldAlert,
  Sparkles,
  Network
} from "lucide-react";
import { useSmartChemAssistant } from "@/hooks/useSmartChemAssistant";

const Assistant = () => {
  // Lift state up to control chat from sidebar
  const { messages, isLoading, sendMessage } = useSmartChemAssistant();

  const quickActions = [
    {
      icon: ShieldAlert,
      label: "ADMET Interpretation",
      prompt: "My molecule flagged as high risk for hERG and DILI toxicity. What structural features usually cause these, and how can I fix them?",
    },
    {
      icon: Network,
      label: "Platform Architecture",
      prompt: "Can you explain how the Hybrid GNN-CNN encoder in SmartChem's architecture fuses sequences and graphs?",
    },
    {
      icon: Lightbulb,
      label: "Medicinal Chemistry",
      prompt: "How can I lower the LogP of my lead compound without losing its aromatic rings or drug-likeness?",
    },
    {
      icon: BrainCircuit,
      label: "Toxicity Filters",
      prompt: "What are PAINS and Brenk alerts in drug discovery, and why do they immediately disqualify a molecule?",
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
                  className="w-full justify-start text-left h-auto py-3 px-4 hover:border-primary/50 transition-colors"
                  onClick={() => sendMessage(action.prompt)}
                  disabled={isLoading}
                >
                  <action.icon className="w-4 h-4 mr-3 text-primary flex-shrink-0" />
                  <span className="text-sm">{action.label}</span>
                </Button>
              </motion.div>
            ))}

            {/* Tips Card */}
            <Card className="p-4 mt-6 bg-secondary/20 border-primary/20">
              <h4 className="font-display text-sm font-semibold text-foreground mb-3 flex items-center gap-2">
                <BrainCircuit className="w-4 h-4 text-primary" />
                Llama-3 Superpowers
              </h4>
              <ul className="space-y-3 text-xs text-muted-foreground">
                <li className="flex items-start gap-2 leading-relaxed">
                  <span className="text-primary mt-0.5">•</span>
                  Ask how to fix specific ADMET violations (like poor BBB penetration)
                </li>
                <li className="flex items-start gap-2 leading-relaxed">
                  <span className="text-primary mt-0.5">•</span>
                  Query about SmartChem's Hybrid VAE and GNN machine learning concepts
                </li>
                <li className="flex items-start gap-2 leading-relaxed">
                  <span className="text-primary mt-0.5">•</span>
                  Learn the chemical rationale behind QED and SAS scoring functions
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
              <ChatInterface
                messages={messages}
                isLoading={isLoading}
                onSendMessage={sendMessage}
              />
            </Card>
          </motion.div>
        </div>
      </div>
    </div>
  );
};

export default Assistant;
