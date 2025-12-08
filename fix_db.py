import asyncio
from backend.database import db

async def run():
    print("Deleting all projects...")
    await db.projects.delete_many({})
    print("Deleting all molecules...")
    await db.molecules.delete_many({})
    print("Done.")

asyncio.run(run())
