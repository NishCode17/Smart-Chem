import * as React from "react";

import { cn } from "@/lib/utils";

const Input = React.forwardRef<HTMLInputElement, React.ComponentProps<"input">>(
  ({ className, type, ...props }, ref) => {
    return (
      <input
        type={type}
        className={cn(
          "flex h-11 w-full rounded-lg border border-primary/20 bg-card/50 backdrop-blur-sm px-4 py-2 text-base text-foreground placeholder:text-muted-foreground transition-all duration-300",
          "focus:outline-none focus:border-primary/50 focus:ring-2 focus:ring-primary/20 focus:shadow-[0_0_15px_hsl(var(--primary)/0.2)]",
          "hover:border-primary/30",
          "file:border-0 file:bg-transparent file:text-sm file:font-medium",
          "disabled:cursor-not-allowed disabled:opacity-50",
          className,
        )}
        ref={ref}
        {...props}
      />
    );
  },
);
Input.displayName = "Input";

export { Input };
