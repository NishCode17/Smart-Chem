import { useState } from 'react'
import axios from 'axios'
import { FlaskConical, Atom, Activity, Zap } from 'lucide-react'

function App() {
  const [numMols, setNumMols] = useState(5)
  const [molecules, setMolecules] = useState([])
  const [loading, setLoading] = useState(false)

  const generateMolecules = async () => {
    setLoading(true)
    try {
      // Connect to your local Python Backend
      const res = await axios.post('http://127.0.0.1:8000/generate', {
        num_molecules: numMols
      })
      setMolecules(res.data.data)
    } catch (err) {
      console.error(err)
      alert("❌ Error: Could not connect to Backend. Is uvicorn running?")
    }
    setLoading(false)
  }

  return (
    <div className="min-h-screen text-slate-800 font-sans">
      {/* 1. TOP NAVIGATION BAR */}
      <nav className="bg-white border-b border-slate-200 px-8 py-4 flex items-center gap-3 sticky top-0 z-50">
        <div className="bg-blue-600 p-2 rounded-lg text-white shadow-lg shadow-blue-200">
          <FlaskConical size={24} strokeWidth={2.5} />
        </div>
        <div>
          <h1 className="text-xl font-bold tracking-tight text-slate-900">SmartChem AI</h1>
          <p className="text-[10px] uppercase tracking-wider text-slate-500 font-bold">De Novo Drug Design Platform</p>
        </div>
        <div className="ml-auto flex gap-4 text-sm font-medium text-slate-500">
          <span className="flex items-center gap-1"><div className="w-2 h-2 rounded-full bg-green-500"></div> System Online</span>
          <span className="flex items-center gap-1"><div className="w-2 h-2 rounded-full bg-blue-500"></div> GPU Active</span>
        </div>
      </nav>

      <main className="p-8 max-w-7xl mx-auto">

        {/* 2. CONTROL PANEL */}
        <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 mb-8 flex flex-col md:flex-row items-center justify-between gap-6">
          <div className="flex-1">
            <h2 className="text-lg font-bold text-slate-900 mb-1 flex items-center gap-2">
              <Activity size={18} className="text-blue-500" /> Generator Configuration
            </h2>
            <p className="text-slate-500 text-sm">Configure the VAE latent space sampling parameters.</p>
          </div>

          <div className="flex items-center gap-6">
            <div className="flex items-center gap-3 bg-slate-50 px-4 py-2 rounded-lg border border-slate-100">
              <label className="text-xs font-bold uppercase text-slate-400">Batch Size</label>
              <input
                type="range" min="1" max="12"
                value={numMols}
                onChange={(e) => setNumMols(parseInt(e.target.value))}
                className="w-32 accent-blue-600 cursor-pointer"
              />
              <span className="font-mono font-bold text-blue-600 text-lg w-8 text-center">{numMols}</span>
            </div>

            <button
              onClick={generateMolecules}
              disabled={loading}
              className="bg-blue-600 hover:bg-blue-700 active:scale-95 transition-all text-white px-8 py-3 rounded-lg font-bold shadow-lg shadow-blue-200 flex items-center gap-2 disabled:opacity-70 disabled:cursor-not-allowed disabled:transform-none"
            >
              {loading ? (
                <>
                  <div className="animate-spin rounded-full h-5 w-5 border-2 border-white border-t-transparent"></div>
                  Synthesizing...
                </>
              ) : (
                <>
                  <Zap size={20} fill="currentColor" />
                  Generate Candidates
                </>
              )}
            </button>
          </div>
        </div>

        {/* 3. RESULTS GRID */}
        {molecules.length > 0 ? (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {molecules.map((mol, idx) => (
              <div key={idx} className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden hover:shadow-xl hover:-translate-y-1 transition-all duration-300 group">

                {/* Header */}
                <div className="px-5 py-3 border-b border-slate-100 flex justify-between items-center bg-slate-50/50">
                  <span className="text-[10px] font-mono font-bold text-slate-400">ID: GEN-{idx + 101}</span>
                  <span className={`text-[10px] uppercase font-bold px-2 py-1 rounded-full border ${mol.properties.status.includes("High")
                      ? "bg-green-50 text-green-700 border-green-200"
                      : "bg-amber-50 text-amber-700 border-amber-200"
                    }`}>
                    {mol.properties.status}
                  </span>
                </div>

                {/* Molecule Image */}
                <div className="h-56 p-4 flex items-center justify-center bg-white relative group-hover:bg-slate-50/30 transition-colors">
                  {mol.properties.image ? (
                    <img src={mol.properties.image} alt="Structure" className="max-h-full max-w-full object-contain mix-blend-multiply filter contrast-125" />
                  ) : (
                    <div className="text-slate-300 flex flex-col items-center">
                      <Atom size={48} />
                      <span className="text-xs mt-2">No 2D Render</span>
                    </div>
                  )}
                </div>

                {/* Metrics */}
                <div className="grid grid-cols-2 border-t border-slate-100 divide-x divide-slate-100">
                  <div className="p-4 flex flex-col items-center">
                    <span className="text-[10px] uppercase tracking-wider text-slate-400 font-bold mb-1">Drug Likeness</span>
                    <span className={`text-xl font-bold ${mol.properties.qed > 0.6 ? 'text-green-600' : 'text-slate-700'}`}>
                      {mol.properties.qed}
                    </span>
                  </div>
                  <div className="p-4 flex flex-col items-center">
                    <span className="text-[10px] uppercase tracking-wider text-slate-400 font-bold mb-1">Solubility (LogP)</span>
                    <span className="text-xl font-bold text-slate-700">{mol.properties.logp}</span>
                  </div>
                </div>

                {/* Footer */}
                <div className="bg-slate-50 px-5 py-3 border-t border-slate-100">
                  <p className="text-[10px] font-mono text-slate-400 truncate" title={mol.properties.smiles}>
                    {mol.properties.smiles}
                  </p>
                </div>
              </div>
            ))}
          </div>
        ) : (
          /* 4. EMPTY STATE */
          !loading && (
            <div className="text-center py-24 border-2 border-dashed border-slate-200 rounded-2xl bg-slate-50/50">
              <div className="bg-white p-4 rounded-full shadow-sm inline-block mb-4">
                <Atom size={48} className="text-blue-500" />
              </div>
              <h3 className="text-lg font-bold text-slate-900">Ready to Synthesize</h3>
              <p className="text-slate-500 max-w-sm mx-auto mt-2">
                The generative model is loaded and ready. Adjust parameters above and click generate to begin discovery.
              </p>
            </div>
          )
        )}
      </main>
    </div>
  )
}

export default App