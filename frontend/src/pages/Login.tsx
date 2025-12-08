import { useState } from "react";
import { motion } from "framer-motion";
import { useNavigate } from "react-router-dom";
import { MoleculeBackground } from "@/components/MoleculeBackground";
import { CyclingText } from "@/components/CyclingText";
import { TransitionAnimation } from "@/components/TransitionAnimation";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { FlaskConical, Lock, User, Eye, EyeOff } from "lucide-react";
import { toast } from "sonner";
import { api } from "@/lib/api";

const Login = () => {
  const navigate = useNavigate();
  const [username, setUsername] = useState("");
  const [email, setEmail] = useState("");
  const [password, setPassword] = useState("");
  const [showPassword, setShowPassword] = useState(false);
  const [isLoading, setIsLoading] = useState(false);
  const [showTransition, setShowTransition] = useState(false);
  const [isRegistering, setIsRegistering] = useState(false);

  const handleAuth = async (e: React.FormEvent) => {
    e.preventDefault();

    if (!username || !password || (isRegistering && !email)) {
      toast.error("Please fill in all credentials");
      return;
    }

    setIsLoading(true);

    try {
      if (isRegistering) {
        await api.register({ username, email, password });
        toast.success("Account created! Logging you in...");
        // Auto-login after register
        await api.login(username, password);
      } else {
        await api.login(username, password);
      }

      localStorage.setItem("smartchem_user", username);
      toast.success(`Welcome back, ${username}!`);

      setShowTransition(true);
    } catch (error: any) {
      toast.error(error.message || "Authentication failed");
      setIsLoading(false);
    }
  };

  const handleTransitionComplete = () => {
    setShowTransition(false);
    navigate("/dashboard");
  };

  return (
    <>
      <TransitionAnimation isAnimating={showTransition} onComplete={handleTransitionComplete} />

      <div className="min-h-screen flex">
        {/* Left Side - Hero Section (70%) */}
        {/* ... existing hero code ... */}
        <div className="hidden lg:flex lg:w-[70%] relative bg-gradient-to-br from-background via-background to-primary/5 overflow-hidden">
          <MoleculeBackground />

          <div className="relative z-10 flex flex-col justify-center px-16 xl:px-24">
            <motion.div
              initial={{ opacity: 0, x: -50 }}
              animate={{ opacity: 1, x: 0 }}
              transition={{ duration: 0.8 }}
            >
              {/* Logo */}
              <div className="flex items-center gap-3 mb-12">
                <div className="w-14 h-14 rounded-2xl bg-gradient-to-br from-primary to-secondary flex items-center justify-center neon-glow">
                  <FlaskConical className="w-8 h-8 text-primary-foreground" />
                </div>
                <span className="font-display text-3xl font-bold text-foreground">
                  SmartChem <span className="gradient-text">AI</span>
                </span>
              </div>

              {/* Cycling text */}
              <CyclingText />

              {/* Tagline */}
              <motion.p
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.4, duration: 0.6 }}
                className="text-xl text-muted-foreground max-w-lg mt-6 leading-relaxed"
              >
                Accelerate drug discovery with AI-powered molecular design,
                analysis, and optimization. From concept to candidate, faster.
              </motion.p>

              {/* Feature highlights */}
              <motion.div
                initial={{ opacity: 0, y: 20 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ delay: 0.6, duration: 0.6 }}
                className="flex gap-6 mt-12"
              >
                {["AI Generation", "3D Analysis", "Lead Optimization"].map((feature, i) => (
                  <div key={i} className="flex items-center gap-2">
                    <div className="w-2 h-2 rounded-full bg-gradient-to-r from-primary to-secondary" />
                    <span className="text-sm text-muted-foreground">{feature}</span>
                  </div>
                ))}
              </motion.div>
            </motion.div>
          </div>

          {/* Decorative elements */}
          <div className="absolute bottom-0 left-0 right-0 h-32 bg-gradient-to-t from-background to-transparent" />
          <div className="absolute top-20 right-20 w-64 h-64 rounded-full bg-primary/10 blur-[100px]" />
          <div className="absolute bottom-40 left-40 w-48 h-48 rounded-full bg-secondary/10 blur-[80px]" />
        </div>

        {/* Right Side - Login Form (30%) */}
        <div className="w-full lg:w-[30%] flex items-center justify-center p-8 bg-background relative">
          {/* Mobile background */}
          <div className="absolute inset-0 lg:hidden overflow-hidden">
            <MoleculeBackground />
            <div className="absolute inset-0 bg-background/80 backdrop-blur-xl" />
          </div>

          <motion.div
            initial={{ opacity: 0, y: 20 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.6, delay: 0.2 }}
            className="w-full max-w-sm relative z-10"
          >
            {/* Mobile logo */}
            <div className="flex items-center justify-center gap-2 mb-8 lg:hidden">
              <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-primary to-secondary flex items-center justify-center">
                <FlaskConical className="w-5 h-5 text-primary-foreground" />
              </div>
              <span className="font-display text-xl font-bold text-foreground">
                SmartChem <span className="gradient-text">AI</span>
              </span>
            </div>

            <Card className="neon-border">
              <CardHeader className="space-y-1 text-center">
                <CardTitle className="font-display text-2xl">
                  {isRegistering ? "Create Account" : "Welcome Back"}
                </CardTitle>
                <p className="text-sm text-muted-foreground">
                  {isRegistering ? "Start your research journey" : "Sign in to your research account"}
                </p>
              </CardHeader>
              <CardContent>
                <form onSubmit={handleAuth} className="space-y-4">
                  <div className="space-y-2">
                    <div className="relative">
                      <User className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-muted-foreground" />
                      <Input
                        type="text"
                        placeholder="Username"
                        value={username}
                        onChange={(e) => setUsername(e.target.value)}
                        className="pl-10"
                      />
                    </div>
                  </div>

                  {isRegistering && (
                    <div className="space-y-2">
                      <div className="relative">
                        <User className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-muted-foreground" />
                        <Input
                          type="email"
                          placeholder="Email Address"
                          value={email}
                          onChange={(e) => setEmail(e.target.value)}
                          className="pl-10"
                        />
                      </div>
                    </div>
                  )}

                  <div className="space-y-2">
                    <div className="relative">
                      <Lock className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-muted-foreground" />
                      <Input
                        type={showPassword ? "text" : "password"}
                        placeholder="Password"
                        value={password}
                        onChange={(e) => setPassword(e.target.value)}
                        className="pl-10 pr-10"
                      />
                      <button
                        type="button"
                        onClick={() => setShowPassword(!showPassword)}
                        className="absolute right-3 top-1/2 -translate-y-1/2 text-muted-foreground hover:text-foreground transition-colors"
                      >
                        {showPassword ? (
                          <EyeOff className="w-4 h-4" />
                        ) : (
                          <Eye className="w-4 h-4" />
                        )}
                      </button>
                    </div>
                  </div>

                  <div className="flex items-center justify-between text-sm">
                    <label className="flex items-center gap-2 text-muted-foreground">
                      <input type="checkbox" className="rounded border-border" />
                      Remember me
                    </label>
                    {!isRegistering && (
                      <a href="#" className="text-primary hover:underline">
                        Forgot password?
                      </a>
                    )}
                  </div>

                  <Button
                    type="submit"
                    variant="neon"
                    size="lg"
                    className="w-full"
                    disabled={isLoading}
                  >
                    {isLoading ? (
                      <motion.div
                        animate={{ rotate: 360 }}
                        transition={{ duration: 1, repeat: Infinity, ease: "linear" }}
                        className="w-5 h-5 border-2 border-primary-foreground border-t-transparent rounded-full"
                      />
                    ) : (
                      isRegistering ? "Register" : "Sign In"
                    )}
                  </Button>
                </form>

                <div className="mt-6 text-center text-sm text-muted-foreground">
                  {isRegistering ? "Already have an account? " : "Don't have an account? "}
                  <button
                    onClick={() => setIsRegistering(!isRegistering)}
                    className="text-primary hover:underline focus:outline-none"
                  >
                    {isRegistering ? "Sign In" : "Register Now"}
                  </button>
                </div>
              </CardContent>
            </Card>

            <p className="text-center text-xs text-muted-foreground mt-6">
              By using SmartChem, you agree to our Terms of Service
            </p>
          </motion.div>
        </div>
      </div>
    </>
  );
};

export default Login;
