import { motion } from "framer-motion";
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  RadarChart,
  PolarGrid,
  PolarAngleAxis,
  PolarRadiusAxis,
  Radar,
} from "recharts";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { AlertTriangle, CheckCircle, XCircle } from "lucide-react";

interface PropertyPanelProps {
  properties: {
    logP: number;
    qed: number;
    molWeight: number;
    hbdCount: number;
    hbaCount: number;
    tpsa: number;
    rotBonds: number;
  };
  admet?: {
    absorption: number;
    distribution: number;
    metabolism: number;
    excretion: number;
    toxicity: number;
  };
  lipinski?: {
    mw: boolean;
    logP: boolean;
    hbd: boolean;
    hba: boolean;
    violations: number;
  };
  activeTab?: string;
}

export const PropertyPanel = ({ properties, admet, lipinski, activeTab = "properties" }: PropertyPanelProps) => {
  const radarData = admet
    ? [
      { property: "Absorption", value: admet.absorption * 100 },
      { property: "Distribution", value: admet.distribution * 100 },
      { property: "Metabolism", value: admet.metabolism * 100 },
      { property: "Excretion", value: admet.excretion * 100 },
      { property: "Toxicity", value: (1 - admet.toxicity) * 100 },
    ]
    : [
      { property: "Absorption", value: 75 },
      { property: "Distribution", value: 68 },
      { property: "Metabolism", value: 82 },
      { property: "Excretion", value: 71 },
      { property: "Toxicity", value: 45 },
    ];

  const barData = [
    { name: "LogP", value: properties.logP, max: 5 },
    { name: "QED", value: properties.qed * 10, max: 10 },
    { name: "TPSA", value: properties.tpsa / 14, max: 10 },
    { name: "RotBonds", value: properties.rotBonds, max: 10 },
  ];

  const lipinskiData = lipinski || {
    mw: properties.molWeight <= 500,
    logP: properties.logP <= 5,
    hbd: properties.hbdCount <= 5,
    hba: properties.hbaCount <= 10,
    violations: 0,
  };

  return (
    <div className="space-y-4">
      {/* Properties Tab Content */}
      {(activeTab === "properties") && (
        <>
          {/* Lipinski Rule of 5 */}
          <motion.div
            initial={{ opacity: 0, x: 20 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ delay: 0.1 }}
          >
            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="text-sm flex items-center gap-2">
                  <CheckCircle className="w-4 h-4 text-primary" />
                  Lipinski Rule of 5
                </CardTitle>
              </CardHeader>
              <CardContent>
                <div className="grid grid-cols-2 gap-2">
                  {[
                    { label: "MW ≤ 500", pass: lipinskiData.mw, value: properties.molWeight.toFixed(0) },
                    { label: "LogP ≤ 5", pass: lipinskiData.logP, value: properties.logP.toFixed(2) },
                    { label: "HBD ≤ 5", pass: lipinskiData.hbd, value: properties.hbdCount },
                    { label: "HBA ≤ 10", pass: lipinskiData.hba, value: properties.hbaCount },
                  ].map((rule, i) => (
                    <div
                      key={i}
                      className={`flex items-center justify-between p-2 rounded-lg text-xs ${rule.pass ? "bg-primary/10" : "bg-destructive/10"
                        }`}
                    >
                      <span className="text-muted-foreground">{rule.label}</span>
                      <span className="flex items-center gap-1">
                        {rule.value}
                        {rule.pass ? (
                          <CheckCircle className="w-3 h-3 text-primary" />
                        ) : (
                          <XCircle className="w-3 h-3 text-destructive" />
                        )}
                      </span>
                    </div>
                  ))}
                </div>
              </CardContent>
            </Card>
          </motion.div>

          {/* Properties Bar Chart */}
          <motion.div
            initial={{ opacity: 0, x: 20 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ delay: 0.3 }}
          >
            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="text-sm">Key Properties</CardTitle>
              </CardHeader>
              <CardContent>
                <ResponsiveContainer width="100%" height={150}>
                  <BarChart data={barData} layout="vertical">
                    <CartesianGrid strokeDasharray="3 3" stroke="hsl(var(--border))" />
                    <XAxis type="number" domain={[0, 10]} tick={{ fill: "hsl(var(--muted-foreground))", fontSize: 10 }} />
                    <YAxis type="category" dataKey="name" tick={{ fill: "hsl(var(--muted-foreground))", fontSize: 10 }} width={60} />
                    <Tooltip
                      contentStyle={{
                        backgroundColor: "hsl(var(--card))",
                        border: "1px solid hsl(var(--border))",
                        borderRadius: "8px",
                      }}
                    />
                    <Bar dataKey="value" fill="url(#gradient)" radius={[0, 4, 4, 0]} />
                    <defs>
                      <linearGradient id="gradient" x1="0" y1="0" x2="1" y2="0">
                        <stop offset="0%" stopColor="hsl(var(--primary))" />
                        <stop offset="100%" stopColor="hsl(var(--secondary))" />
                      </linearGradient>
                    </defs>
                  </BarChart>
                </ResponsiveContainer>
              </CardContent>
            </Card>
          </motion.div>
        </>
      )}

      {/* ADMET Tab Content */}
      {activeTab === "admet" && (
        <motion.div
          initial={{ opacity: 0, x: 20 }}
          animate={{ opacity: 1, x: 0 }}
          transition={{ delay: 0.2 }}
        >
          <Card>
            <CardHeader className="pb-2">
              <CardTitle className="text-sm">ADMET Profile</CardTitle>
            </CardHeader>
            <CardContent>
              <ResponsiveContainer width="100%" height={300}>
                <RadarChart data={radarData}>
                  <PolarGrid stroke="hsl(var(--border))" />
                  <PolarAngleAxis
                    dataKey="property"
                    tick={{ fill: "hsl(var(--muted-foreground))", fontSize: 11 }}
                  />
                  <PolarRadiusAxis
                    angle={30}
                    domain={[0, 100]}
                    tick={{ fill: "hsl(var(--muted-foreground))", fontSize: 10 }}
                  />
                  <Radar
                    name="ADMET"
                    dataKey="value"
                    stroke="hsl(var(--primary))"
                    fill="hsl(var(--primary))"
                    fillOpacity={0.3}
                  />
                </RadarChart>
              </ResponsiveContainer>
              <div className="mt-4 space-y-2">
                <div className="text-xs text-muted-foreground bg-secondary/20 p-3 rounded-lg">
                  High absorption and good metabolic stability predicted. Low toxicity risk indicated by AI model confidence of 89%.
                </div>
              </div>
            </CardContent>
          </Card>
        </motion.div>
      )}

      {/* Toxicity Tab Content */}
      {activeTab === "toxicity" && (
        <motion.div
          initial={{ opacity: 0, x: 20 }}
          animate={{ opacity: 1, x: 0 }}
          transition={{ delay: 0.4 }}
        >
          <Card>
            <CardHeader className="pb-2">
              <CardTitle className="text-sm flex items-center gap-2">
                <AlertTriangle className="w-4 h-4 text-yellow-500" />
                Toxicity Alerts
              </CardTitle>
            </CardHeader>
            <CardContent>
              <div className="space-y-2">
                {[
                  { label: "PAINS Filter", status: "pass" },
                  { label: "Brenk Alerts", status: "warning" },
                  { label: "NIH Alerts", status: "pass" },
                  { label: "Liver Toxicity", status: "pass" }, // Added extra sample data
                  { label: "Mutagenicity", status: "pass" },
                ].map((alert, i) => (
                  <div key={i} className="flex items-center justify-between text-sm py-2 border-b border-border/50 last:border-0">
                    <span className="text-muted-foreground">{alert.label}</span>
                    <span
                      className={`px-2 py-0.5 rounded-full text-xs font-medium ${alert.status === "pass"
                          ? "bg-green-500/10 text-green-500"
                          : "bg-yellow-500/10 text-yellow-500"
                        }`}
                    >
                      {alert.status === "pass" ? "Pass" : "Review"}
                    </span>
                  </div>
                ))}
              </div>
              <div className="mt-4 p-3 bg-yellow-500/5 border border-yellow-500/20 rounded-lg">
                <p className="text-xs text-yellow-600">
                  <strong>Note:</strong> One structural alert detected (Brenk). Specific substructure may require optimization to reduce potential reactivity.
                </p>
              </div>
            </CardContent>
          </Card>
        </motion.div>
      )}
    </div>
  );
};
