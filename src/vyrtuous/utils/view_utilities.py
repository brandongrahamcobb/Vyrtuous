def limit_available_to_top_25_by_member_count(available):
    available.sort(key=lambda a: getattr(a, "member_count", 0), reverse=True)
    top_25 = available[:25]
    top_ids = {a.id for a in top_25}
    filtered = []
    seen = set()
    for a_list in available.values():
        for a in a_list:
            if a.id in top_ids and a.id not in seen:
                seen.add(a.id)
                filtered.append(a)
    return filtered
