import { useState, useRef, useEffect } from "react";
import { motion, AnimatePresence } from "framer-motion";
import { Send, Bot, User, Sparkles, Loader2 } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { useSmartChemAssistant } from "@/hooks/useSmartChemAssistant";

interface Message {
  id: string;
  role: "user" | "assistant";
  content: string;
  timestamp: Date;
}

interface ChatInterfaceProps {
  moleculeContext?: {
    smiles: string;
    properties: Record<string, number>;
  };
  compact?: boolean;
}

export const ChatInterface = ({ moleculeContext, compact = false }: ChatInterfaceProps) => {
  /* Hook Integration */
  const { messages, isLoading, sendMessage } = useSmartChemAssistant();
  const [input, setInput] = useState("");
  const messagesEndRef = useRef<HTMLDivElement>(null);

  const scrollToBottom = () => {
    messagesEndRef.current?.scrollIntoView({ behavior: "smooth" });
  };

  useEffect(() => {
    scrollToBottom();
  }, [messages]);

  const handleSend = async () => {
    if (!input.trim() || isLoading) return;

    const currentInput = input;
    setInput("");

    // Pass context if available
    await sendMessage(currentInput, moleculeContext);
  };

  return (
    <div className={`flex flex-col h-full ${compact ? "max-h-[400px]" : ""}`}>
      {/* Header */}
      {/* Header Removed for cleaner layout */}

      {/* Messages */}
      <div className="flex-1 overflow-y-auto p-4 space-y-4 scrollbar-thin">
        <AnimatePresence>
          {messages.map((message) => (
            <motion.div
              key={message.id}
              initial={{ opacity: 0, y: 10 }}
              animate={{ opacity: 1, y: 0 }}
              exit={{ opacity: 0, y: -10 }}
              className={`flex gap-3 ${message.role === "user" ? "flex-row-reverse" : ""}`}
            >
              <div
                className={`w-8 h-8 rounded-full flex items-center justify-center flex-shrink-0 ${message.role === "assistant"
                  ? "bg-gradient-to-br from-primary to-secondary"
                  : "bg-muted"
                  }`}
              >
                {message.role === "assistant" ? (
                  <Bot className="w-4 h-4 text-primary-foreground" />
                ) : (
                  <User className="w-4 h-4 text-foreground" />
                )}
              </div>
              <div
                className={`max-w-[80%] p-3 rounded-xl ${message.role === "assistant" ? "chat-bubble-ai" : "chat-bubble-user"
                  }`}
              >
                <p className="text-sm text-foreground leading-relaxed">{message.content}</p>
                <span className="text-xs text-muted-foreground mt-2 block">
                  {message.timestamp.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" })}
                </span>
              </div>
            </motion.div>
          ))}
        </AnimatePresence>

        {isLoading && (
          <motion.div
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            className="flex gap-3"
          >
            <div className="w-8 h-8 rounded-full bg-gradient-to-br from-primary to-secondary flex items-center justify-center">
              <Bot className="w-4 h-4 text-primary-foreground" />
            </div>
            <div className="chat-bubble-ai p-3 rounded-xl">
              <Loader2 className="w-4 h-4 animate-spin text-primary" />
            </div>
          </motion.div>
        )}
        <div ref={messagesEndRef} />
      </div>

      {/* Context indicator */}
      {moleculeContext && (
        <div className="px-4 py-2 border-t border-border bg-primary/5">
          <p className="text-xs text-muted-foreground">
            <Sparkles className="w-3 h-3 inline mr-1" />
            Context: Analyzing molecule {moleculeContext.smiles.slice(0, 20)}...
          </p>
        </div>
      )}

      {/* Input */}
      <div className="p-4 border-t border-border">
        <div className="flex gap-2">
          <Input
            value={input}
            onChange={(e) => setInput(e.target.value)}
            onKeyPress={(e) => e.key === "Enter" && handleSend()}
            placeholder="Ask Dr. SmartChem..."
            className="flex-1"
          />
          <Button
            variant="neon"
            size="icon"
            onClick={handleSend}
            disabled={isLoading || !input.trim()}
          >
            <Send className="w-4 h-4" />
          </Button>
        </div>
      </div>
    </div>
  );
};
