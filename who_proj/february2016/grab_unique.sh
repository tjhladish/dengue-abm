while read tag; do
  grep -l $tag auto_output/newjunk.err* | head -n 1
done <alt_tags
