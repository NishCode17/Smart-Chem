
import os
import json
import requests
from pydantic import BaseModel
from typing import Optional, Dict, Any
from dotenv import load_dotenv

load_dotenv()

GROQ_API_KEY = os.environ.get("GROQ_API_KEY")
GROQ_API_URL = "https://api.groq.com/openai/v1/chat/completions"
MODEL_NAME = "llama-3.3-70b-versatile"

class ChatRequest(BaseModel):
    message: str
    context: Optional[Dict[str, Any]] = None

def build_context_string(context: Optional[Dict[str, Any]]) -> str:
    if not context:
        return "No specific molecule context provided."
    
    lines = ["Current Molecule Data:"]
    
    for key, value in context.items():
        formatted_key = key.replace("_", " ").title()
        if isinstance(value, dict):
            lines.append(f"- {formatted_key}:")
            for sub_k, sub_v in value.items():
                lines.append(f"  * {sub_k.replace('_', ' ').title()}: {sub_v}")
        elif isinstance(value, list):
            lines.append(f"- {formatted_key}: {', '.join(map(str, value))}")
        elif isinstance(value, (str, int, float, bool)):
            lines.append(f"- {formatted_key}: {value}")
            
    return "\n".join(lines)

def get_groq_response(message: str, context: Optional[Dict[str, Any]]) -> str:
    if not GROQ_API_KEY:
        return "Error: GROQ_API_KEY environment variable not set."

    context_str = build_context_string(context)
    
    system_prompt = (
        "You are Dr. SmartChem, an AI medicinal chemistry assistant. "
        "Your role is to explain molecular properties, analyze ADMET behavior, assess Lipinski compliance, "
        "identify toxicity alerts, provide structural reasoning, and suggest optimization ideas.\n\n"
        "RULES:\n"
        "1. NEVER hallucinate structural details. Base your reasoning ONLY on the provided context.\n"
        "2. If the context is missing specific data, state that you don't have that information.\n"
        "3. Be professional, precise, and helpful to a medicinal chemist.\n"
        "4. Keep responses concise but informative."
    )

    user_content = f"Context:\n{context_str}\n\nUser Question:\n{message}"

    payload = {
        "model": MODEL_NAME,
        "messages": [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_content}
        ],
        "temperature": 0.5,
        "max_tokens": 1024
    }

    try:
        response = requests.post(
            GROQ_API_URL,
            headers={
                "Authorization": f"Bearer {GROQ_API_KEY}",
                "Content-Type": "application/json"
            },
            json=payload,
            timeout=30
        )
        response.raise_for_status()
        data = response.json()
        return data["choices"][0]["message"]["content"]
    except Exception as e:
        print(f"Groq API Error: {e}")
        return "I apologize, but I'm having trouble connecting to my knowledge base right now. Please try again later."
