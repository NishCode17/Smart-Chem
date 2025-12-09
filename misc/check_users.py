
import asyncio
import os
from motor.motor_asyncio import AsyncIOMotorClient
from dotenv import load_dotenv

load_dotenv()

from passlib.context import CryptContext

pwd_context = CryptContext(schemes=["argon2"], deprecated="auto")

async def check_users():
    # ... (existing connection code) ...
    mongo_url = os.getenv("MONGO_URI")
    db_name = os.getenv("MONGO_DB")
    client = AsyncIOMotorClient(mongo_url)
    db = client[db_name]

    # Test Hashing
    test_pass = "password123"
    hashed = pwd_context.hash(test_pass)
    print(f"Test Hash: {hashed}")
    verified = pwd_context.verify(test_pass, hashed)
    print(f"Test Verify: {verified}")
    
    # ... (connection above) ...

    # fetch users
    try:
        cursor = db.users.find().sort("created_at", -1).limit(5)
        users = await cursor.to_list(length=5)
        
        print(f"Found {len(users)} users.")
        for user in users:
            print(f"Username: {user.get('username')}")
            print(f"Email: {user.get('email')}")
            # print first 20 chars of hash to verify it exists and looks like a hash
            p_hash = user.get('password_hash', 'NO_HASH_FOUND')
            print(f"Hash: {p_hash}")
            if p_hash.startswith("$argon2"):
                print("Type: Argon2")
            elif p_hash.startswith("$2b$"):
                 print("Type: Bcrypt")
            else:
                 print("Type: Unknown")
            print("-" * 20)
            
            # Try verifying 'password' (just a guess for debugging)
            try:
                if pwd_context.verify("password", p_hash):
                    print("  [DEBUG] Password is 'password'")
                if pwd_context.verify("123456", p_hash):
                    print("  [DEBUG] Password is '123456'")
            except Exception as e:
                print(f"  [DEBUG] Verify Error: {e}")

    except Exception as e:
        print(f"Error fetching users: {e}")

if __name__ == "__main__":
    asyncio.run(check_users())
