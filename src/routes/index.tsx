import { createFileRoute } from "@tanstack/react-router";
import { WindLab } from "@/components/wind-lab";

export const Route = createFileRoute("/")({ component: Home });

function Home() {
  return (
    <main className="min-h-screen bg-bg text-fg">
      <WindLab />
    </main>
  );
}
