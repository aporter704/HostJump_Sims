#accessing sqlite3 in python
import sqlite3

#create working directory "hostjumps.db"
con = sqlite3.connect("hostjumps.db")

#create database cursor to fetch **TABLE NAME**
cur = con.cursor() 

#create database table to store **TABLE NAME**, called "**TABLE NAME**"
#columns include X, Y, Z

cur.execute("CREATE TABLE **TABLE NAME**(tree_name, tmrca, jumps, first_migration_age)")

#verify table has been made, should return "true" if created, "none" if something went wrong

res = cur.execute("SELECT name FROM sqlite_master")
res.fetchone()

#example of importing **TABLE NAME** into table via each row (eg. name for each tree and corresponding **TABLE NAME** for tmrca, and number of jumps, for one row will be labelled "tree1", 2019.1, 3)
cur.execute("""INSERT INTO **TABLE NAME** VALUES ("tree1", 2019.1, 3), ("tree2", 2018.5, 2)""")

#commit changes to table
con.commit()

#verify **TABLE NAME** have been loaded into table
res = cur.execute("SELECT tmrca FROM **TABLE NAME**")

#importing many values into table example - need to use the correct number of ? placeholders for each column
data =[(###**TABLE NAME** from trees in comma delimited form)]
cur.executemany("INSERT INTO **TABLE NAME** VALUES(?,?,?)", data)
con.commit()

#can use SELECT to pick out, order, or otherwise interact with table
res = new_cur.execute("SELECT tmrca, jumps, tree_name FROM **TABLE NAME** ORDER BY tmrca DESC"):
#selects all rows
res.fetchall()
#selects many rows
res.fetchmany(size=x)
#select one row, with certain parameters
tree_name, tmrca,jumps = res.fetchone()
print(f'The oldest tmrca of trees is {tree_name!r}, with a tmrca of {tmrca}, and {jumps} host-jumps')


#to insert more data as a script
con.executescript("""
INSERT INTO **TABLE NAME** (tree_name, tmrca, jumps) VALUES ('tree5', '2013', '5')
INSERT INTO **TABLE NAME** (tree_name, tmrca, jumps) VALUES ('tree6', '2017', '1')
INSERT INTO **TABLE NAME** (tree_name, tmrca, jumps) VALUES ('tree7', '2016', '0')
""")

#to close connection 
close()