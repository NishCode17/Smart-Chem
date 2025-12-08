import { Navigate, useLocation } from "react-router-dom";
import { api } from "@/lib/api";

interface RequireAuthProps {
    children: JSX.Element;
}

const RequireAuth: React.FC<RequireAuthProps> = ({ children }) => {
    const auth = api.isAuthenticated(); // We need to add this method to api.ts first
    const location = useLocation();

    if (!auth) {
        return <Navigate to="/" state={{ from: location }} replace />;
    }

    return children;
};

export default RequireAuth;
