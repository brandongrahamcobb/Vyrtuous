import requests
from typing import List, Dict, Any

class RobinhoodCryptoTracker:
    BASE_URL = "https://api.robinhood.com/"

    def __init__(self, auth_token: str):
        self.auth_token = auth_token
        self.headers = {
            "Authorization": f"Bearer {self.auth_token}",
            "Content-Type": "application/json",
        }

    def get_crypto_list(self) -> List[Dict[str, Any]]:
        url = f"{self.BASE_URL}cryptocurrencies/"
        response = requests.get(url, headers=self.headers)
        response.raise_for_status()
        return response.json().get("results", [])

    def get_crypto_quote(self, symbol: str) -> Dict[str, Any]:
        url = f"{self.BASE_URL}marketdata/forex/quotes/{symbol}/"
        response = requests.get(url, headers=self.headers)
        response.raise_for_status()
        return response.json()

    def track_statistics(self, symbols: List[str]) -> Dict[str, Dict[str, Any]]:
        stats = {}
        for symbol in symbols:
            try:
                quote = self.get_crypto_quote(symbol)
                stats[symbol] = {
                    "price": quote.get("mark_price"),
                    "bid": quote.get("bid_price"),
                    "ask": quote.get("ask_price"),
                    "volume": quote.get("volume"),
                }
            except requests.exceptions.RequestException as e:
                stats[symbol] = {"error": str(e)}
        return stats

if __name__ == "__main__":
    tracker = RobinhoodCryptoTracker(auth_token="your_auth_token_here")
    crypto_symbols = ["BTC", "ETH", "DOGE"]
    statistics = tracker.track_statistics(crypto_symbols)
    for symbol, data in statistics.items():
        print(f"{symbol}: {data}")
