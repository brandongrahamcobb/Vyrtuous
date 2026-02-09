def limit_available_to_top_25_by_member_count(available):
    all_key = "all"
    items = []
    if all_key in available:
        items.append(all_key)
    objects = []
    for k, v in available.items():
        if k == all_key:
            continue
        if isinstance(v, (list, tuple, set)):
            objects.extend(v)
        else:
            objects.append(v)
    objects.sort(key=lambda a: getattr(a, "member_count", 0), reverse=True)
    top_25 = objects[:25]
    return items + top_25
