import { motion } from "framer-motion";
import { Button } from "@/components/ui/button";
import { Card, CardContent } from "@/components/ui/card";
import { FlaskConical, Save, Send, Sparkles } from "lucide-react";

interface MoleculeCardProps {
  id: string;
  image?: string;
  smiles: string;
  properties: {
    logP: number;
    qed: number;
    molWeight: number;
    hbdCount: number;
    hbaCount: number;
  };
  onSave?: () => void;
  onSendToLab?: () => void;
}

export const MoleculeCard = ({
  id,
  image,
  smiles,
  properties,
  onSave,
  onSendToLab,
}: MoleculeCardProps) => {
  return (
    <motion.div
      initial={{ opacity: 0, scale: 0.95 }}
      animate={{ opacity: 1, scale: 1 }}
      whileHover={{ y: -4 }}
      transition={{ duration: 0.3 }}
    >
      <Card className="overflow-hidden hover-glow">
        <div className="h-48 bg-white flex items-center justify-center relative overflow-hidden">
          {image ? (
            <img src={image} alt="Molecule Structure" className="w-full h-full object-contain p-2" />
          ) : (
            <div className="absolute inset-0 flex items-center justify-center">
              <div className="w-20 h-20 rounded-full bg-secondary/30 flex items-center justify-center">
                <FlaskConical className="w-8 h-8 text-primary opacity-50" />
              </div>
            </div>
          )}
        </div>

        <CardContent className="p-4 space-y-4">
          {/* SMILES string */}
          <div className="font-mono text-xs text-muted-foreground truncate bg-muted/50 p-2 rounded-lg">
            {smiles}
          </div>

          {/* Properties grid */}
          <div className="grid grid-cols-2 gap-2 text-sm">
            <div className="flex justify-between">
              <span className="text-muted-foreground">LogP</span>
              <span className="font-medium text-foreground">{properties.logP.toFixed(2)}</span>
            </div>
            <div className="flex justify-between">
              <span className="text-muted-foreground">QED</span>
              <span className="font-medium text-foreground">{properties.qed.toFixed(2)}</span>
            </div>
            <div className="flex justify-between">
              <span className="text-muted-foreground">MW</span>
              <span className="font-medium text-foreground">{properties.molWeight.toFixed(1)}</span>
            </div>
            <div className="flex justify-between">
              <span className="text-muted-foreground">HBD/HBA</span>
              <span className="font-medium text-foreground">
                {properties.hbdCount}/{properties.hbaCount}
              </span>
            </div>
          </div>

          {/* Action buttons */}
          <div className="flex gap-2 relative z-10">
            <Button
              variant="outline"
              size="sm"
              className="flex-1 cursor-pointer hover:bg-secondary"
              onClick={(e) => {
                e.stopPropagation();
                onSave?.();
              }}
            >
              <Save className="w-4 h-4 mr-1" />
              Save
            </Button>
            <Button
              variant="neon"
              size="sm"
              className="flex-1 cursor-pointer"
              onClick={(e) => {
                e.stopPropagation();
                onSendToLab?.();
              }}
            >
              <Send className="w-4 h-4 mr-1" />
              To Lab
            </Button>
          </div>
        </CardContent>
      </Card>
    </motion.div>
  );
};
