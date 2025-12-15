import sys

def check_asyncpg():
    try:
        import asyncpg
        print(f"✅ asyncpg is installed: {asyncpg.__version__}")
        return True
    except ImportError as e:
        print(f"❌ asyncpg not found: {e}")
        return False

if __name__ == "__main__":
    success = check_asyncpg()
    sys.exit(0 if success else 1)
