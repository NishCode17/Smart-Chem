
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
    history: Optional[list] = []

def build_context_string(context: Optional[Dict[str, Any]]) -> Optional[str]:
    if not context:
        return None
    
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

# Load readme context
try:
    readme_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "README.md")
    with open(readme_path, "r", encoding="utf-8") as f:
        PROJECT_KNOWLEDGE_BASE = f.read()
except Exception as e:
    PROJECT_KNOWLEDGE_BASE = "Project knowledge base could not be loaded."

def get_groq_response(message: str, context: Optional[Dict[str, Any]], history: Optional[list] = None) -> str:
    if not GROQ_API_KEY:
        return "Error: GROQ_API_KEY environment variable not set."

    context_str = build_context_string(context)
    
    system_prompt = (
        "You are Dr. SmartChem, the highly advanced AI medicinal chemistry assistant for the 'SmartChem' de novo drug discovery platform. "
        "You were created as a B.Tech Major Project by students Nishant Singh Bisht, Navya Singh, Nishant Kumar Singh, and Shivam Sahani "
        "from the Department of Computer Science and Engineering at United College of Engineering & Research (UCER), Prayagraj.\n\n"
        
        "=== SMARTCHEM PROJECT KNOWLEDGE BASE ===\n"
        f"{PROJECT_KNOWLEDGE_BASE}\n"
        "==========================================\n\n"
        
        "YOUR BEHAVIOR & CAPABILITIES:\n"
        "1. You have complete knowledge of the SmartChem platform based on the Knowledge Base above. If asked about the platform's architecture, results, or features, use that information confidently.\n"
        "2. You possess vast pre-trained knowledge of medicinal chemistry literature, mechanisms of action, and drug design principles.\n"
        "3. CRITICAL: If asked to 'search literature' or find papers, DO NOT say 'I cannot search the internet' or 'I don't have real-time access'. Instead, act as an expert and RECALL foundational papers, well-known studies, and established literature from your training data. Cite real authors and concepts.\n"
        "4. If a 'Context (Molecule Data)' block is provided by the user, explain its molecular properties, analyze ADMET behavior, assess Lipinski compliance, and suggest optimization ideas. NEVER hallucinate specific molecular properties—rely ONLY on the provided Context block.\n"
        "5. If NO Context block is provided, simply engage in a natural conversation as an expert chemist. You can freely answer questions, recommend papers, and discuss chemistry topics. DO NOT complain about missing context unless the user specifically asks you to analyze a molecule.\n"
        "6. Be professional, highly technical, and extremely helpful to a medicinal chemist or a university Viva examiner."
    )

    if context_str:
        user_content = f"Context (Molecule Data):\n{context_str}\n\nUser Question:\n{message}"
    else:
        user_content = message

    messages_payload = [{"role": "system", "content": system_prompt}]
    
    if history:
        for msg in history[-10:]:  # Keep last 10 messages
            if isinstance(msg, dict) and "role" in msg and "content" in msg:
                messages_payload.append({"role": msg["role"], "content": msg["content"]})
                
    messages_payload.append({"role": "user", "content": user_content})

    print("=== DEBUG PAYLOAD ===")
    for idx, m in enumerate(messages_payload):
        print(f"[{idx}] {m['role'].upper()}: {m['content'][:100]}...")
    print("=====================")

    payload = {
        "model": MODEL_NAME,
        "messages": messages_payload,

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
