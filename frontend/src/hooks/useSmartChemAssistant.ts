
import { useState } from 'react';

export interface Message {
    id: string;
    role: "user" | "assistant";
    content: string;
    timestamp: Date;
}

interface MoleculeContext {
    smiles: string;
    properties: Record<string, number | string | boolean>;
}

export const useSmartChemAssistant = () => {
    const [messages, setMessages] = useState<Message[]>([
        {
            id: "1",
            role: "assistant",
            content: "Hello! I'm Dr. SmartChem, your AI chemistry assistant. How can I help you today? I can analyze molecules, suggest optimizations, or answer questions about drug design.",
            timestamp: new Date(),
        },
    ]);
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);

    const sendMessage = async (content: string, context?: MoleculeContext) => {
        if (!content.trim()) return;

        const userMessage: Message = {
            id: Date.now().toString(),
            role: "user",
            content,
            timestamp: new Date(),
        };

        setMessages((prev) => [...prev, userMessage]);
        setIsLoading(true);
        setError(null);

        try {
            const response = await fetch("http://localhost:8000/assistant/chat", {
                method: "POST",
                headers: {
                    "Content-Type": "application/json",
                },
                body: JSON.stringify({
                    message: content,
                    context: context ? {
                        smiles: context.smiles,
                        ...context.properties
                    } : undefined,
                }),
            });

            if (!response.ok) {
                throw new Error("Failed to get response from assistant");
            }

            const data = await response.json();

            const assistantMessage: Message = {
                id: (Date.now() + 1).toString(),
                role: "assistant",
                content: data.reply,
                timestamp: new Date(),
            };

            setMessages((prev) => [...prev, assistantMessage]);
        } catch (err) {
            console.error(err);
            setError("Failed to connect to Dr. SmartChem. Please try again.");
            // Optionally add an error message to the chat
            const errorMessage: Message = {
                id: (Date.now() + 1).toString(),
                role: "assistant",
                content: "I'm having trouble connecting to the server. Please check your connection and try again.",
                timestamp: new Date(),
            };
            setMessages((prev) => [...prev, errorMessage]);
        } finally {
            setIsLoading(false);
        }
    };

    return {
        messages,
        isLoading,
        error,
        sendMessage,
    };
};
