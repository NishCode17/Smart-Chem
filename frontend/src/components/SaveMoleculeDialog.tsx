import React, { useEffect, useState } from "react";

import {
    Dialog,
    DialogContent,
    DialogHeader,
    DialogTitle,
    DialogDescription,
} from "@/components/ui/dialog";

import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";

import {
    Select,
    SelectTrigger,
    SelectContent,
    SelectItem,
    SelectValue,
} from "@/components/ui/select";

import { api, Project_ } from "@/lib/api";
import { toast } from "sonner";

// ------------------------------------------------------------
// NAMED EXPORT ONLY
// ------------------------------------------------------------
export function SaveMoleculeDialog({
    open,
    onOpenChange,
    molecule,
    defaultName = "New Molecule",
    sourceType = "generator",
}) {
    const [projects, setProjects] = useState<Project_[]>([]);
    const [selectedProject, setSelectedProject] = useState<string>("");
    const [name, setName] = useState(defaultName);
    const [isSaving, setIsSaving] = useState(false);

    // Load projects when dialog opens
    useEffect(() => {
        if (open) {
            setName(defaultName);
            loadProjects();
        }
    }, [open, defaultName]);

    const loadProjects = async () => {
        try {
            const p = await api.getProjects();
            setProjects(p);

            // Fix controlled/uncontrolled issue
            if (p.length > 0) {
                setSelectedProject((prev) => prev || p[0]._id);
            }
        } catch (err) {
            console.error(err);
            toast.error("Failed to load projects");
        }
    };

    const handleSave = async () => {
        if (!selectedProject) {
            toast.error("Please select a project");
            return;
        }

        if (!molecule?.smiles) {
            toast.error("No molecule to save");
            return;
        }

        setIsSaving(true);

        try {
            await api.saveMolecule({
                name,
                project_id: selectedProject,
                smiles: molecule.smiles,
                generated_by: sourceType,
                tags: ["Lead"],

                properties: {
                    logp: molecule.properties?.logP ?? 0,
                    qed: molecule.properties?.qed ?? 0,
                    mw: molecule.properties?.molWeight ?? 0,
                    tpsa: molecule.properties?.tpsa ?? 0,
                    hbd: molecule.properties?.hbdCount ?? 0,
                    hba: molecule.properties?.hbaCount ?? 0,
                    rot_bonds: molecule.properties?.rotatable ?? 0,
                },

                // VALIDATION FIX: Ensure ADMET object has correct shape
                // If molecule.admet has 'mw' (descriptors), use default scores instead
                admet: (molecule.admet && 'absorption' in molecule.admet) ? molecule.admet : {
                    absorption: 0,
                    distribution: 0,
                    metabolism: 0,
                    excretion: 0,
                    toxicity: 0,
                },

                tox_alerts: molecule.tox_alerts ?? {},
            });

            toast.success("Molecule saved!");
            onOpenChange(false);
        } catch (err: any) {
            console.error(err);
            toast.error(err.message || "Failed to save molecule");
        } finally {
            setIsSaving(false);
        }
    };

    return (
        <Dialog open={open} onOpenChange={onOpenChange}>
            <DialogContent>
                <DialogHeader>
                    <DialogTitle>Save Molecule</DialogTitle>
                    <DialogDescription>
                        Store this molecule inside a project folder.
                    </DialogDescription>
                </DialogHeader>

                {/* Molecule name */}
                <div className="space-y-2">
                    <label className="text-sm font-medium">Name</label>
                    <Input value={name} onChange={(e) => setName(e.target.value)} />
                </div>

                {/* Project select */}
                <div className="space-y-2 mt-4">
                    <label className="text-sm font-medium">Save To Project</label>

                    <Select
                        value={selectedProject}
                        onValueChange={setSelectedProject}
                    >
                        <SelectTrigger>
                            <SelectValue placeholder="Choose a project" />
                        </SelectTrigger>

                        <SelectContent>
                            {projects.map((p) => (
                                <SelectItem key={p._id} value={p._id}>
                                    {p.name}
                                </SelectItem>
                            ))}
                        </SelectContent>
                    </Select>
                </div>

                <Button
                    className="w-full mt-6"
                    disabled={isSaving}
                    onClick={handleSave}
                >
                    {isSaving ? "Saving..." : "Save Molecule"}
                </Button>
            </DialogContent>
        </Dialog>
    );
}
