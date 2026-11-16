import { createRoot } from "react-dom/client";
import { WindLab } from "@/components/wind-lab";
import "@/styles.css";

const root = document.getElementById("root");
if (root) {
  createRoot(root).render(<WindLab />);
}
