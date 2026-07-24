import yaml

from vyrtuous.models import permission


def test_load_groups():
    groups = permission.load_groups("/app/vyrtuous/groups.yml")
    permissions = permission.resolve_group_permissions("Administrator", groups)
    assert "command.listing.blacklists" in permissions
