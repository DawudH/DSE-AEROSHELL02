n=open("library.bib", "r")
bib=n.readlines()
n.close()

nr=3

for i in range(len(bib)):
    if  bib[i].find('author =')==0:
        l=bib[i]
        authors=l.split('and')
        if len(authors)>nr:
            l=''
            for j in range(nr):
                l= l +authors[j] +' and '
                print l
            l=l +'{ Others}} \n,'
            print l
            bib[i]=l
            

f=open("Bibliography.bib", 'w')
for i in range(len(bib)):
    f.write(bib[i])
f.close

