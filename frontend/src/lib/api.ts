const API_BASE_URL = "http://localhost:8000";

// --- Types ---
export interface AuthResponse {
    access_token: string;
    refresh_token: string;
    token_type: string;
}

export interface UserCreate {
    username: string;
    email: string;
    password: string;
}

export interface Project_ {
    _id: string;                        // Only Mongo _id allowed
    name: string;
    description?: string;
    created_at: string;
}

export interface Molecule {
    id: string;
    project_id: string;
    name: string;
    smiles: string;
    generated_by: string;
    tags: string[];
    properties?: {
        logp: number;
        qed: number;
        mw: number;
        tpsa: number;
        hbd: number;
        hba: number;
        rot_bonds: number;
    };
    admet?: {
        absorption: number;
        distribution: number;
        metabolism: number;
        excretion: number;
        toxicity: number;
    };
    tox_alerts?: {
        pains: boolean;
        brenk: boolean;
        nih: boolean;
    };
    created_at: string;
}

export interface GenerateRequest {
    num_molecules: number;
    target_qed?: number;
    target_logp?: number;
    target_sas?: number;
}

export interface OptimizeRequest {
    smiles: string;
    target_qed: number;
    target_logp: number;
}

export interface MoleculeData {
    properties: {
        smiles: string;
        image?: string;
        logp: number;
        qed: number;
        admet?: {
            mw: number;
            hbd: number;
            hba: number;
            tpsa: number;
        };
        [key: string]: any;
    };
}

// --- Auth Utils ---
function getAccessToken(): string | null {
    return localStorage.getItem("access_token");
}

function getRefreshToken(): string | null {
    return localStorage.getItem("refresh_token");
}

function setTokens(access: string, refresh: string) {
    localStorage.setItem("access_token", access);
    localStorage.setItem("refresh_token", refresh);
}

function clearTokens() {
    localStorage.removeItem("access_token");
    localStorage.removeItem("refresh_token");
    localStorage.removeItem("smartchem_user");
}

// --- Client ---
async function fetchWithAuth<T>(endpoint: string, options: RequestInit = {}): Promise<T> {
    let token = getAccessToken();

    const headers: HeadersInit = {
        "Content-Type": "application/json",
        ...options.headers,
    };

    if (token) {
        headers["Authorization"] = `Bearer ${token}`;
    }

    let res = await fetch(`${API_BASE_URL}${endpoint}`, {
        ...options,
        headers,
    });

    // Handle 401 — Refresh Token
    if (res.status === 401) {
        const refreshToken = getRefreshToken();
        if (refreshToken) {
            try {
                const refreshRes = await fetch(`${API_BASE_URL}/auth/refresh`, {
                    method: "POST",
                    headers: { "Content-Type": "application/json" },
                    body: JSON.stringify({ refresh_token: refreshToken }),
                });

                if (refreshRes.ok) {
                    const data: AuthResponse = await refreshRes.json();
                    setTokens(data.access_token, data.refresh_token);

                    headers["Authorization"] = `Bearer ${data.access_token}`;
                    res = await fetch(`${API_BASE_URL}${endpoint}`, {
                        ...options,
                        headers,
                    });
                } else {
                    clearTokens();
                    window.location.href = "/";
                    throw new Error("Session expired");
                }
            } catch (err) {
                clearTokens();
                window.location.href = "/";
                throw err;
            }
        } else {
            clearTokens();
        }
    }

    if (!res.ok) {
        const errorData = await res.json().catch(() => ({ detail: res.statusText }));
        let errorMessage = errorData.detail || `API error: ${res.status}`;
        if (typeof errorMessage === "object") {
            errorMessage = JSON.stringify(errorMessage, null, 2);
        }
        throw new Error(errorMessage);
    }

    return res.json();
}

// --- API Methods ---
export const api = {
    // Auth
    login: async (username: string, password: string): Promise<AuthResponse> => {
        const formData = new URLSearchParams();
        formData.append("username", username);
        formData.append("password", password);

        const res = await fetch(`${API_BASE_URL}/auth/login`, {
            method: "POST",
            headers: { "Content-Type": "application/x-www-form-urlencoded" },
            body: formData,
        });

        if (!res.ok) {
            const err = await res.json();
            throw new Error(err.detail || "Login failed");
        }

        const data: AuthResponse = await res.json();
        setTokens(data.access_token, data.refresh_token);
        localStorage.setItem("smartchem_user", username);

        return data;
    },

    register: async (user: UserCreate) =>
        fetch(`${API_BASE_URL}/auth/register`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(user),
        }).then(async (res) => {
            if (!res.ok) {
                const err = await res.json();
                throw new Error(err.detail || "Registration failed");
            }
            return res.json();
        }),

    logout: () => {
        clearTokens();
        window.location.href = "/";
    },

    // Projects
    getProjects: () =>
        fetchWithAuth<any[]>("/projects/").then((projects) =>
            projects.map((p) => ({
                _id: p.id || p._id,         // FIXED: Backend returns 'id'
                name: p.name,
                description: p.description,
                created_at: p.created_at,
            }))
        ),

    createProject: (name: string, description?: string) =>
        fetchWithAuth<any>("/projects/", {
            method: "POST",
            body: JSON.stringify({ name, description }),
        }).then((p) => ({
            _id: p.id || p._id,          // FIXED: Backend returns 'id'
            name: p.name,
            description: p.description,
            created_at: p.created_at,
        })),

    deleteProject: (id: string) =>
        fetchWithAuth(`/projects/${id}`, { method: "DELETE" }),

    // Molecules
    getMolecules: (projectId?: string) => {
        const query = projectId ? `?project_id=${projectId}` : "";
        return fetchWithAuth<Molecule[]>(`/molecules/${query}`);
    },

    getMolecule: (id: string) =>
        fetchWithAuth<Molecule>(`/molecules/${id}`),

    saveMolecule: (mol: Partial<Molecule>) =>
        fetchWithAuth<Molecule>("/molecules/", {
            method: "POST",
            body: JSON.stringify({
                ...mol,
                project_id: mol.project_id,   // ENSURED CORRECT FIELD
                generated_by: mol.generated_by || "imported",
            }),
        }),

    deleteMolecule: (id: string) =>
        fetchWithAuth(`/molecules/${id}`, { method: "DELETE" }),

    // AI / Utils
    generate: (req: GenerateRequest) =>
        fetchWithAuth<{ data: any[] }>("/generate", {
            method: "POST",
            body: JSON.stringify(req),
        }),

    generateTargeted: (req: GenerateRequest) =>
        fetchWithAuth<{ data: any[] }>("/generate/targeted", {
            method: "POST",
            body: JSON.stringify(req),
        }),

    optimizeLead: (req: OptimizeRequest) =>
        fetchWithAuth<{ data: any[] }>("/optimize/lead", {
            method: "POST",
            body: JSON.stringify(req),
        }),

    get3DStructure: (smiles: string) =>
        fetchWithAuth<{ mol_block: string }>("/utils/3d", {
            method: "POST",
            body: JSON.stringify({ smiles }),
        }),

    isAuthenticated: () => {
        return !!getAccessToken();
    },
};
