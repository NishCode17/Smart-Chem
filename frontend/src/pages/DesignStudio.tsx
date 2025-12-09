import { useState } from "react";
import { motion } from "framer-motion";
import { PageHeader } from "@/components/PageHeader";
import { MoleculeCard } from "@/components/MoleculeCard";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import {
  Wand2,
  Shuffle,
  Target,
  Loader2,
  Sparkles,
  Settings2
} from "lucide-react";
import { toast } from "sonner";
import { useNavigate } from "react-router-dom";
import { api, MoleculeData } from "@/lib/api";

interface GeneratedMolecule {
  id: string;
  image?: string;
  smiles: string;
  properties: {
    logP: number;
    qed: number;
    molWeight: number;
    hbdCount: number;
    hbaCount: number;
    tpsa: number;
    rotatable: number;
  };
  admet?: any;
  tox_alerts?: any;
}

import { SaveMoleculeDialog } from "@/components/SaveMoleculeDialog";

const DesignStudio = () => {
  // ... (lines 37-51 in original are fine) ...
  const navigate = useNavigate();
  const [activeTab, setActiveTab] = useState("generator");
  const [isGenerating, setIsGenerating] = useState(false);

  // State for inputs
  const [inputSmiles, setInputSmiles] = useState("");
  const [targetQED, setTargetQED] = useState(0.8);
  const [targetLogP, setTargetLogP] = useState(2.5);

  const [generatedMolecules, setGeneratedMolecules] = useState<GeneratedMolecule[]>([]);

  // Save Dialog State
  const [saveDialogMolecule, setSaveDialogMolecule] = useState<GeneratedMolecule | null>(null);
  const [isSaveDialogOpen, setIsSaveDialogOpen] = useState(false);

  // Helper to map API data to UI format
  const mapResults = (results: MoleculeData[]) => {
    return results.map((r, i) => ({
      id: `${Date.now()}-${i}`,
      image: r.properties.image,
      smiles: r.properties.smiles,
      properties: {
        logP: r.properties.logp,
        qed: r.properties.qed,
        molWeight: (r.properties as any).admet?.mw || 0,
        hbdCount: (r.properties as any).admet?.hbd || 0,
        hbaCount: (r.properties as any).admet?.hba || 0,
        tpsa: (r.properties as any).admet?.tpsa || 0,
        rotatable: (r.properties as any).admet?.rotatable || 0,
      },
      admet: (r.properties as any).admet,
      tox_alerts: (r.properties as any).tox_alerts
    }));
  };

  const handleRandomGenerate = async () => {
    setIsGenerating(true);
    try {
      const res = await api.generate({ num_molecules: 3 });
      setGeneratedMolecules(mapResults(res.data));
      toast.success(`Generated ${res.data.length} new molecules`);
    } catch (e: any) {
      toast.error(e.message || "Generation failed");
    } finally {
      setIsGenerating(false);
    }
  };

  const handleTargetedGenerate = async () => {
    setIsGenerating(true);
    try {
      // Use Property Optimization from scratch
      const res = await api.generateTargeted({
        num_molecules: 3,
        target_qed: targetQED,
        target_logp: targetLogP,
        target_sas: 3.0 // Default internal
      });
      setGeneratedMolecules(mapResults(res.data));
      toast.success("Generated targeted molecules");
    } catch (e: any) {
      toast.error(e.message || "Generation failed");
    } finally {
      setIsGenerating(false);
    }
  };

  const handleOptimize = async () => {
    if (!inputSmiles) {
      toast.error("Please enter a SMILES string to optimize");
      return;
    }
    setIsGenerating(true);
    try {
      const res = await api.optimizeLead({
        smiles: inputSmiles,
        target_qed: targetQED,
        target_logp: targetLogP
      });
      setGeneratedMolecules(mapResults(res.data));
      toast.success("Lead optimization complete!");
    } catch (e: any) {
      toast.error(e.message || "Optimization failed");
    } finally {
      setIsGenerating(false);
    }
  };

  const handleSave = (molecule: GeneratedMolecule) => {
    setSaveDialogMolecule(molecule);
    setIsSaveDialogOpen(true);
  };

  const handleSendToLab = (molecule: GeneratedMolecule) => {
    // Map to the format VirtualLab/API expects (lowercase property keys)
    const labMolecule = {
      id: molecule.id,
      smiles: molecule.smiles,
      name: `Generated-${molecule.id.slice(-4)}`,
      properties: {
        logp: molecule.properties.logP, // Fix case mismatch
        qed: molecule.properties.qed,
        mw: molecule.properties.molWeight,
        tpsa: molecule.properties.tpsa,
        hbd: molecule.properties.hbdCount,
        hba: molecule.properties.hbaCount,
        rot_bonds: molecule.properties.rotatable
      },
      admet: molecule.admet,
      tox_alerts: molecule.tox_alerts,
      generated_by: activeTab === "optimizer" ? "optimizer" : "generator"
    };

    navigate("/virtual-lab", { state: { molecule: labMolecule } });
  };

  return (
    <div className="min-h-screen bg-background p-6">
      <SaveMoleculeDialog
        open={isSaveDialogOpen}
        onOpenChange={setIsSaveDialogOpen}
        molecule={saveDialogMolecule}
        defaultName={saveDialogMolecule ? `Compound-${saveDialogMolecule.id.slice(-4)}` : "New Molecule"}
        sourceType={activeTab === "optimizer" ? "optimizer" : "generator"}
      />

      <div className="container mx-auto max-w-7xl">
        <PageHeader
          title="Design Studio"
          description="Generate and optimize novel drug candidates"
        />

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Left Panel - Controls */}
          <div className="lg:col-span-1 space-y-4">
            <Tabs value={activeTab} onValueChange={setActiveTab}>
              <TabsList className="w-full grid grid-cols-2">
                <TabsTrigger value="generator">Generator</TabsTrigger>
                <TabsTrigger value="optimizer">Optimizer</TabsTrigger>
              </TabsList>

              <TabsContent value="generator" className="space-y-4 mt-4">
                <Card>
                  <CardHeader className="pb-3">
                    <CardTitle className="text-sm flex items-center gap-2">
                      <Shuffle className="w-4 h-4 text-primary" />
                      Random Generation
                    </CardTitle>
                  </CardHeader>
                  <CardContent>
                    <p className="text-xs text-muted-foreground mb-4">
                      Generate novel drug-like molecules from scratch using AI
                    </p>
                    <Button
                      variant="neon"
                      className="w-full"
                      onClick={handleRandomGenerate}
                      disabled={isGenerating}
                    >
                      {isGenerating ? (
                        <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                      ) : (
                        <Sparkles className="w-4 h-4 mr-2" />
                      )}
                      Generate Random
                    </Button>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader className="pb-3">
                    <CardTitle className="text-sm flex items-center gap-2">
                      <Target className="w-4 h-4 text-primary" />
                      Targeted Properties
                    </CardTitle>
                  </CardHeader>
                  <CardContent className="space-y-3">
                    <p className="text-xs text-muted-foreground">
                      Generate molecules with specific properties
                    </p>
                    <div className="space-y-2">
                      <div className="flex justify-between text-xs">
                        <label>Target QED</label>
                        <span className="text-muted-foreground">{targetQED}</span>
                      </div>
                      <Input
                        type="number"
                        step="0.1"
                        min="0"
                        max="1"
                        value={targetQED}
                        onChange={(e) => setTargetQED(parseFloat(e.target.value))}
                        className="font-mono text-xs"
                      />
                      <div className="flex justify-between text-xs">
                        <label>Target LogP</label>
                        <span className="text-muted-foreground">{targetLogP}</span>
                      </div>
                      <Input
                        type="number"
                        step="0.5"
                        value={targetLogP}
                        onChange={(e) => setTargetLogP(parseFloat(e.target.value))}
                        className="font-mono text-xs"
                      />
                    </div>
                    <Button
                      variant="outline"
                      className="w-full"
                      onClick={handleTargetedGenerate}
                      disabled={isGenerating}
                    >
                      {isGenerating ? (
                        <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                      ) : (
                        <Target className="w-4 h-4 mr-2" />
                      )}
                      Generate Targeted
                    </Button>
                  </CardContent>
                </Card>
              </TabsContent>

              <TabsContent value="optimizer" className="space-y-4 mt-4">
                <Card>
                  <CardHeader className="pb-3">
                    <CardTitle className="text-sm flex items-center gap-2">
                      <Wand2 className="w-4 h-4 text-primary" />
                      Lead Optimizer
                    </CardTitle>
                  </CardHeader>
                  <CardContent className="space-y-3">
                    <p className="text-xs text-muted-foreground">
                      Optimize a lead compound for better drug-likeness and activity
                    </p>
                    <Input
                      placeholder="Enter lead SMILES..."
                      value={inputSmiles}
                      onChange={(e) => setInputSmiles(e.target.value)}
                      className="font-mono text-xs"
                    />
                    <div className="flex items-center gap-2 text-xs text-muted-foreground">
                      <Settings2 className="w-3 h-3" />
                      Using Magic Wand AI
                    </div>
                    <Button
                      variant="neon"
                      className="w-full"
                      onClick={handleOptimize}
                      disabled={isGenerating}
                    >
                      {isGenerating ? (
                        <Loader2 className="w-4 h-4 mr-2 animate-spin" />
                      ) : (
                        <Wand2 className="w-4 h-4 mr-2" />
                      )}
                      Optimize Lead
                    </Button>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader className="pb-3">
                    <CardTitle className="text-sm">Optimization Targets</CardTitle>
                  </CardHeader>
                  <CardContent>
                    <div className="space-y-4">
                      <div className="space-y-1">
                        <label className="text-xs font-medium">Drug-likeness (QED)</label>
                        <Input
                          type="number"
                          step="0.1"
                          value={targetQED}
                          onChange={(e) => setTargetQED(parseFloat(e.target.value))}
                          className="h-8 text-xs font-mono"
                        />
                      </div>
                      <div className="space-y-1">
                        <label className="text-xs font-medium">Lipophilicity (LogP)</label>
                        <Input
                          type="number"
                          step="0.5"
                          value={targetLogP}
                          onChange={(e) => setTargetLogP(parseFloat(e.target.value))}
                          className="h-8 text-xs font-mono"
                        />
                      </div>
                    </div>
                  </CardContent>
                </Card>
              </TabsContent>
            </Tabs>
          </div>

          {/* Right Panel - Results */}
          <div className="lg:col-span-2">
            <motion.div
              initial={{ opacity: 0 }}
              animate={{ opacity: 1 }}
              className="space-y-4"
            >
              <div className="flex items-center justify-between">
                <h2 className="font-display text-lg font-semibold text-foreground">
                  Generated Molecules
                </h2>
                <span className="text-sm text-muted-foreground">
                  {generatedMolecules.length} results
                </span>
              </div>

              {generatedMolecules.length === 0 ? (
                <Card className="p-12 text-center">
                  <Sparkles className="w-12 h-12 mx-auto text-primary/30 mb-4" />
                  <p className="text-muted-foreground">
                    No molecules generated yet. Use the tools on the left to create new candidates.
                  </p>
                </Card>
              ) : (
                <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 gap-4">
                  {generatedMolecules.map((molecule) => (
                    <MoleculeCard
                      key={molecule.id}
                      {...molecule}
                      onSave={() => handleSave(molecule)}
                      onSendToLab={() => handleSendToLab(molecule)}
                    />
                  ))}
                </div>
              )}
            </motion.div>
          </div>
        </div>
      </div>
    </div>
  );
};


export default DesignStudio;
