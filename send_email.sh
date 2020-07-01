#using mailx
#-s subject
#-c CC address
#-r reply-to "name<email>"
#-a attachment

#message body by echo
echo "message body" | mailx -s "test email from rescomp" -r "NGS_pipeline" recipient1,recipient2

#message body by file
mailx -s "test email from rescomp" -r "NGS_pipeline" recipient1,recipient2 < body.txt
cat body.txt | mailx -s "test email from rescomp" -r "NGS_pipeline" recipient1,recipient2
