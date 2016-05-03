#Makefile
main: main.o parameter.o initial.o core_acc.o gas_acc.o typeI_migration.o typeII_migration.o trap.o next.o gimpact.o output.o
	gcc -o a.out main.o parameter.o initial.o core_acc.o gas_acc.o typeI_migration.o typeII_migration.o trap.o next.o gimpact.o output.o

main.o: main.c
	gcc -c main.c

parameter.o: parameter.c
	gcc -c parameter.c

initial.o: initial.c
	gcc -c initial.c

core_acc.o: core_acc.c
	gcc -c core_acc.c

gac_acc.o: gas_acc.c
	gcc -c gas_acc.c

typeI_migration.o: typeI_migration.c
	gcc -c typeI_migration.c

typeII_migration.o: typeII_migration.c
	gcc -c typeII_migration.c

trap.o: trap.c
	gcc -c trap.c

next.o: next.c
	gcc -c next.c

gimpact.o: gimpact.c
	gcc -c gimpact.c

output.o: output.c
	gcc -c output.c

clean:
	rm -f *.o
