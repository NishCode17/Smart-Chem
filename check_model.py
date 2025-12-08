from backend.models import ProjectDB

print("Creating ProjectDB...")

p = ProjectDB(name="Test Project", user_id="123")
print(p)

print("Generated ID =", p.id)
