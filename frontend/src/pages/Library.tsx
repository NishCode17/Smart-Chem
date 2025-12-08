import { useState, useEffect } from "react";
import { motion, AnimatePresence } from "framer-motion";
import { useNavigate } from "react-router-dom";
import { PageHeader } from "@/components/PageHeader";
import { Card, CardContent } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu";
import {
  Dialog,
  DialogContent,
  DialogDescription,
  DialogFooter,
  DialogHeader,
  DialogTitle,
  DialogTrigger,
} from "@/components/ui/dialog";
import {
  Folder,
  FolderPlus,
  Search,
  Grid,
  List,
  MoreVertical,
  FlaskConical,
  Trash2,
  Edit,
  Download,
  Tag,
  ChevronRight,
} from "lucide-react";
import { toast } from "sonner";
import { api, Project_, Molecule } from "@/lib/api";

const LibraryPage = () => {
  const navigate = useNavigate();
  const [viewMode, setViewMode] = useState<"grid" | "list">("grid");
  const [searchQuery, setSearchQuery] = useState("");
  const [selectedProjectId, setSelectedProjectId] = useState<string | null>(null);
  const [showNewProjectDialog, setShowNewProjectDialog] = useState(false);
  const [newProjectName, setNewProjectName] = useState("");
  const [newProjectDesc, setNewProjectDesc] = useState("");

  const [projects, setProjects] = useState<Project_[]>([]);
  const [molecules, setMolecules] = useState<Molecule[]>([]);
  const [isLoadingProjects, setIsLoadingProjects] = useState(false);
  const [isLoadingMolecules, setIsLoadingMolecules] = useState(false);

  // Initial Load
  useEffect(() => {
    fetchProjects();
  }, []);

  // Fetch Molecules when Project Changes
  useEffect(() => {
    if (selectedProjectId) {
      fetchMolecules(selectedProjectId);
    } else {
      setMolecules([]);
    }
  }, [selectedProjectId]);

  const fetchProjects = async () => {
    setIsLoadingProjects(true);
    try {
      const data = await api.getProjects();
      setProjects(data);
      if (data.length > 0 && !selectedProjectId) {
        setSelectedProjectId(data[0]._id); // Auto-select first
      }
    } catch (error) {
      toast.error("Failed to load projects");
    } finally {
      setIsLoadingProjects(false);
    }
  };

  const fetchMolecules = async (projectId: string) => {
    setIsLoadingMolecules(true);
    try {
      const data = await api.getMolecules(projectId);
      setMolecules(data);
    } catch (error) {
      toast.error("Failed to load molecules");
    } finally {
      setIsLoadingMolecules(false);
    }
  };

  const handleCreateProject = async () => {
    if (!newProjectName.trim()) return;
    try {
      const newProject = await api.createProject(newProjectName, newProjectDesc);
      setProjects([...projects, newProject]);
      setSelectedProjectId(newProject._id);
      setNewProjectName("");
      setNewProjectDesc("");
      setShowNewProjectDialog(false);
      toast.success(`Project "${newProject.name}" created`);
    } catch (error) {
      toast.error("Failed to create project");
    }
  };

  const handleDeleteMolecule = async (id: string, e: React.MouseEvent) => {
    e.stopPropagation();
    if (!confirm("Are you sure you want to delete this molecule?")) return;
    try {
      await api.deleteMolecule(id);
      setMolecules(molecules.filter(m => m.id !== id));
      toast.success("Molecule deleted");
    } catch (error) {
      toast.error("Failed to delete molecule");
    }
  };

  const handleDeleteProject = async (id: string, e: React.MouseEvent) => {
    e.stopPropagation();
    if (!confirm("Delete project and all its molecules?")) return;
    try {
      await api.deleteProject(id);
      setProjects(projects.filter(p => p._id !== id));
      if (selectedProjectId === id) setSelectedProjectId(projects[0]?._id || null);
      toast.success("Project deleted");
    } catch (error) {
      toast.error("Failed to delete project");
    }
  };

  const filteredMolecules = molecules.filter(
    (mol) =>
      mol.name.toLowerCase().includes(searchQuery.toLowerCase()) ||
      mol.smiles.toLowerCase().includes(searchQuery.toLowerCase()) ||
      mol.tags?.some((tag) => tag.toLowerCase().includes(searchQuery.toLowerCase()))
  );

  const handleOpenMolecule = (molecule: Molecule) => {
    navigate("/virtual-lab", {
      state: {
        molecule: {
          ...molecule,
          // Map backend simplified properties to frontend expectations if needed
          // Assuming backend structure matches or is sufficient
        },
      },
    });
  };

  const handleExportCSV = () => {
    // Basic CSV Export
    const headers = ["Name", "SMILES", "Tags", "QED", "LogP", "Created At"];
    const rows = molecules.map(m => [
      m.name,
      m.smiles,
      m.tags.join(";"),
      m.properties?.qed || "",
      m.properties?.logp || "",
      new Date(m.created_at).toLocaleDateString()
    ]);

    const csvContent = "data:text/csv;charset=utf-8,"
      + [headers.join(","), ...rows.map(e => e.join(","))].join("\n");

    const encodedUri = encodeURI(csvContent);
    const link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", `molecules_${selectedProjectId}.csv`);
    document.body.appendChild(link);
    link.click();
  };

  return (
    <div className="min-h-screen bg-background p-6">
      <div className="container mx-auto max-w-7xl">
        <PageHeader
          title="Library"
          description="Organize and manage your molecule collection"
        />

        {/* Toolbar */}
        <motion.div
          initial={{ opacity: 0, y: 10 }}
          animate={{ opacity: 1, y: 0 }}
          className="flex flex-col md:flex-row gap-4 mb-6"
        >
          <div className="relative flex-1">
            <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-muted-foreground" />
            <Input
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              placeholder="Search molecules by name, SMILES, or tags..."
              className="pl-10"
            />
          </div>
          <div className="flex gap-2">
            <Dialog open={showNewProjectDialog} onOpenChange={setShowNewProjectDialog}>
              <DialogTrigger asChild>
                <Button variant="outline">
                  <FolderPlus className="w-4 h-4 mr-2" />
                  New Project
                </Button>
              </DialogTrigger>
              <DialogContent>
                <DialogHeader>
                  <DialogTitle>Create New Project</DialogTitle>
                  <DialogDescription>
                    Organize your molecules into project folders
                  </DialogDescription>
                </DialogHeader>
                <div className="space-y-4">
                  <Input
                    value={newProjectName}
                    onChange={(e) => setNewProjectName(e.target.value)}
                    placeholder="Project Name..."
                  />
                  <Input
                    value={newProjectDesc}
                    onChange={(e) => setNewProjectDesc(e.target.value)}
                    placeholder="Description (Optional)"
                  />
                </div>
                <DialogFooter>
                  <Button variant="outline" onClick={() => setShowNewProjectDialog(false)}>
                    Cancel
                  </Button>
                  <Button variant="neon" onClick={handleCreateProject}>
                    Create
                  </Button>
                </DialogFooter>
              </DialogContent>
            </Dialog>
            <Button variant="outline" onClick={handleExportCSV} disabled={!selectedProjectId || molecules.length === 0}>
              <Download className="w-4 h-4 mr-2" />
              Export
            </Button>
            <div className="flex border border-border rounded-lg overflow-hidden">
              <Button
                variant={viewMode === "grid" ? "default" : "ghost"}
                size="icon"
                onClick={() => setViewMode("grid")}
                className="rounded-none"
              >
                <Grid className="w-4 h-4" />
              </Button>
              <Button
                variant={viewMode === "list" ? "default" : "ghost"}
                size="icon"
                onClick={() => setViewMode("list")}
                className="rounded-none"
              >
                <List className="w-4 h-4" />
              </Button>
            </div>
          </div>
        </motion.div>

        <div className="grid grid-cols-1 lg:grid-cols-4 gap-6">
          {/* Projects Sidebar */}
          <motion.div
            initial={{ opacity: 0, x: -20 }}
            animate={{ opacity: 1, x: 0 }}
            className="lg:col-span-1 space-y-2"
          >
            <div className="flex items-center justify-between mb-4">
              <h3 className="font-display text-sm font-semibold text-foreground flex items-center gap-2">
                <Folder className="w-4 h-4 text-primary" />
                Projects
              </h3>
            </div>

            {isLoadingProjects ? (
              <div className="text-sm text-muted-foreground animate-pulse">Loading projects...</div>
            ) : projects.length === 0 ? (
              <div className="text-sm text-muted-foreground">No projects yet. Create one!</div>
            ) : (
              projects.map((project) => (
                <div key={project._id} className="group flex items-center gap-1">
                  <Button
                    variant={selectedProjectId === project._id ? "secondary" : "ghost"}
                    className="flex-1 justify-start truncate"
                    onClick={() => setSelectedProjectId(project._id)}
                  >
                    <Folder className={`w-4 h-4 mr-2 ${selectedProjectId === project._id ? "text-primary" : "text-muted-foreground"}`} />
                    <span className="truncate">{project.name}</span>
                  </Button>
                  <Button variant="ghost" size="icon" className="h-8 w-8 opacity-0 group-hover:opacity-100 transition-opacity text-destructive hover:text-destructive" onClick={(e) => handleDeleteProject(project._id, e)}>
                    <Trash2 className="w-3 h-3" />
                  </Button>
                </div>
              ))
            )}
          </motion.div>

          {/* Molecules Grid/List */}
          <motion.div
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            className="lg:col-span-3"
          >
            {isLoadingMolecules ? (
              <div className="flex items-center justify-center p-12">
                <motion.div
                  animate={{ rotate: 360 }}
                  transition={{ duration: 1, repeat: Infinity, ease: "linear" }}
                  className="w-8 h-8 border-2 border-primary border-t-transparent rounded-full"
                />
              </div>
            ) : filteredMolecules.length === 0 ? (
              <Card className="p-12 text-center border-dashed">
                <FlaskConical className="w-12 h-12 mx-auto text-primary/30 mb-4" />
                <p className="text-muted-foreground">
                  {selectedProjectId ? "No molecules in this project." : "Select a project to view molecules."}
                </p>
              </Card>
            ) : viewMode === "grid" ? (
              <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 gap-4">
                <AnimatePresence>
                  {filteredMolecules.map((molecule, i) => (
                    <motion.div
                      key={molecule.id}
                      initial={{ opacity: 0, scale: 0.95 }}
                      animate={{ opacity: 1, scale: 1 }}
                      exit={{ opacity: 0, scale: 0.95 }}
                      transition={{ delay: i * 0.05 }}
                    >
                      <Card
                        className="cursor-pointer hover-glow transition-all group"
                        onClick={() => handleOpenMolecule(molecule)}
                      >
                        {/* Placeholder visual */}
                        <div className="h-24 bg-gradient-to-br from-muted to-background/50 flex items-center justify-center relative overflow-hidden">
                          <div className="text-[10px] font-mono text-muted-foreground/30 absolute inset-0 p-2 break-all leading-none opacity-50">
                            {molecule.smiles}
                          </div>
                          <FlaskConical className="w-8 h-8 text-primary/20 z-10" />
                        </div>

                        <CardContent className="p-4">
                          <div className="flex items-start justify-between mb-2">
                            <h4 className="font-semibold text-foreground truncate pr-2" title={molecule.name}>{molecule.name}</h4>
                            <DropdownMenu>
                              <DropdownMenuTrigger asChild onClick={(e) => e.stopPropagation()}>
                                <Button variant="ghost" size="icon" className="h-6 w-6 -mt-1 -mr-2">
                                  <MoreVertical className="w-3 h-3" />
                                </Button>
                              </DropdownMenuTrigger>
                              <DropdownMenuContent align="end">
                                <DropdownMenuItem onClick={(e) => handleDeleteMolecule(molecule.id, e)} className="text-destructive">
                                  <Trash2 className="w-4 h-4 mr-2" />
                                  Delete
                                </DropdownMenuItem>
                              </DropdownMenuContent>
                            </DropdownMenu>
                          </div>

                          <div className="flex items-center justify-between mt-3 text-xs">
                            <div className="flex flex-wrap gap-1">
                              {molecule.tags && molecule.tags.map((tag) => (
                                <span key={tag} className="px-1.5 py-0.5 rounded-full bg-primary/10 text-primary border border-primary/20">
                                  {tag}
                                </span>
                              ))}
                            </div>
                          </div>

                          <div className="mt-4 grid grid-cols-3 gap-2 text-xs text-muted-foreground border-t border-border/50 pt-2">
                            {molecule.properties && (
                              <>
                                <div className="text-center">
                                  <div className="font-semibold text-foreground">{molecule.properties.qed.toFixed(2)}</div>
                                  <div className="text-[10px]">QED</div>
                                </div>
                                <div className="text-center border-l border-border/50">
                                  <div className="font-semibold text-foreground">{molecule.properties.logp.toFixed(1)}</div>
                                  <div className="text-[10px]">LogP</div>
                                </div>
                                <div className="text-center border-l border-border/50">
                                  <div className="font-semibold text-foreground">{molecule.properties.mw.toFixed(0)}</div>
                                  <div className="text-[10px]">MW</div>
                                </div>
                              </>
                            )}
                          </div>
                        </CardContent>
                      </Card>
                    </motion.div>
                  ))}
                </AnimatePresence>
              </div>
            ) : (
              <Card className="overflow-hidden">
                <div className="divide-y divide-border">
                  {filteredMolecules.map((molecule) => (
                    <div
                      key={molecule.id}
                      className="flex items-center p-3 hover:bg-muted/50 cursor-pointer transition-colors"
                      onClick={() => handleOpenMolecule(molecule)}
                    >
                      <div className="w-8 h-8 rounded bg-primary/10 flex items-center justify-center mr-3">
                        <FlaskConical className="w-4 h-4 text-primary" />
                      </div>
                      <div className="flex-1 min-w-0">
                        <h4 className="font-medium text-sm text-foreground truncate">{molecule.name}</h4>
                        <p className="font-mono text-xs text-muted-foreground truncate opacity-70">
                          {molecule.smiles}
                        </p>
                      </div>
                      <div className="text-xs text-muted-foreground mx-4">
                        {new Date(molecule.created_at).toLocaleDateString()}
                      </div>
                      <Button size="icon" variant="ghost" className="h-8 w-8 text-destructive/50 hover:text-destructive" onClick={(e) => handleDeleteMolecule(molecule.id, e)}>
                        <Trash2 className="w-4 h-4" />
                      </Button>
                    </div>
                  ))}
                </div>
              </Card>
            )}
          </motion.div>
        </div>
      </div>
    </div>
  );
};

export default LibraryPage;
