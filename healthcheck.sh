HEARTBEAT=$(cat /tmp/vyrtuous_heartbeat)

if [ "$?" == "1" ] || [ "$HEARTBEAT" == "1" ]; then
  exit 1
elif [ "$HEARTBEAT" == "0" ]; then
  exit 0
fi
