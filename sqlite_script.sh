
#accessing sqlite3 in python
import sqlite3

#create working directory "hostjumps.db"
con = sqlite3.connect("hostjumps.db")

#create database cursor to fetch results
cur = con.cursor() 

#create database table to store results, called "results"
#columns include X, Y, Z

cur.execute("CREATE TABLE results(sim_name, tmrca, jumps)")

#verify table has been made, should return "true" if created, "none" if something went wrong

res = cur.execute("SELECT name FROM sqlite_master")
res.fetchone()

#example of importing results into table via each row (eg. name for each simulation and corresponding results for tmrca, and number of jumps, for one row will be labelled "simulation1", 2019.1, 3)
cur.execute("""INSERT INTO results VALUES ("simulation1", 2019.1, 3), ("simulation2", 2018.5, 2)""")

#commit changes to table
con.commit()

#verify results have been loaded into table
res = cur.execute("SELECT tmrca FROM results")

#importing many values into table example - need to use the correct number of ? placeholders for each column
data =[(###results from simulations in comma delimited form)]
cur.executemany("INSERT INTO results VALUES(?,?,?)", data)
con.commit()

#can use SELECT to pick out, order, or otherwise interact with table
res = new_cur.execute("SELECT tmrca, jumps, sim_name FROM results ORDER BY tmrca DESC"):
#selects all rows
res.fetchall()
#selects many rows
res.fetchmany(size=x)
#select one row, with certain parameters
sim_name, tmrca,jumps = res.fetchone()
print(f'The oldest tmrca of simultations is {sim_name!r}, with a tmrca of {tmrca}, and {jumps} host-jumps')


#to insert more data as a script
con.executescript("""
INSERT INTO results (sim_name, tmrca, jumps) VALUES ('simulation5', '2013', '5')
INSERT INTO results (sim_name, tmrca, jumps) VALUES ('simulation6', '2017', '1')
INSERT INTO results (sim_name, tmrca, jumps) VALUES ('simulation7', '2016', '0')
""")

#to close connection 
close()