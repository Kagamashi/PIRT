import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import java.util.HashMap;

public class PirtMaster {

    //To Lepiej zostaw
    public static ArrayList<Connection> connections_;
    public static boolean AWAIT = true;
    public static double X = -1;
    //Do tego momentu



    //jak jest na true to wyswietla ladne macierze, jesli na false wyswietla macierze do skopiowania i wklejenia do programu
    public static boolean LadneMacierze=true;
    public static boolean pokaz_iteracje = true; //pokazuje iteracje erlanga

    enum TypZadania {
        obl_AZ_Z_S, //Oblicza Macierz AZ z S
        obl_AZ_Z_S_z_Roznymi_beta, //Oblicza Macierz AZ z S z roznymi beta
        obl_AZ_Z_Macierzy, // Oblicza AZ Z Macierzy WZ i A

        obl_Zasoby_Wiazek, //Oblicza zasoby wiązek
        obl_Pokaz_Wszystkie_Polaczenia_Przechodzace_Przez_Polaczenie, //Wyswietla wszystkie mozliwe polaczenie przechodzace przez dane polaczenie

        obl_Ilosc_Stanowisk_Prawdopodobienstwo_Oczekiwania, //Oblicza ilosc stanowisk przy podanym minimalnym prawdopodobienstwie oczekiwania
        obl_Erlang_2_dla_N_A, // Oblicza Druga Funkcje Erlanga
        obl_N_dla_progu_Erlang_2, // Oblicza N Dla progu z drugiej Funkcji Erlanga
        obl_Erlang_dla_N_A,// Oblicza Pierwsza Funkcje Erlanga
        obl_N_dla_progu_Erlang,// Oblicza N Dla progu z pierwszej Funkcji Erlanga

    }

    //Typ Zadania
    static TypZadania typ = TypZadania.obl_Ilosc_Stanowisk_Prawdopodobienstwo_Oczekiwania;

    //AZ Z S
    public static float[] S = {1500,2000,0,0,0,0,3000,2500};
    public static float beta = 0.06f;
    public static float[] betaLista = {0.015f,0.02f,0,0,0.025f};

    //AZ Z Macierzy
    public static float[] Azm = {1290,1530,2550};
    public static double[][] WZ_M ={
            {0.3, 0.1, 0.2, X, 0.2, 0.2},
            {0.1, 0.2, 0.1, X, 0.3, 0.3},
            {0.2, 0.2, 0.2, X, 0.2, 0.2},
            {X,     X,   X, X,   X,   X},
            {0.1, 0.1, 0.1, X, 0.3, 0.4},
            {0.3, 0.1, 0.2, X, 0.1, 0.3},
    };

    //obl_Zasoby_Wiazek
    public static int[] polaczenie = {4, 5}; //polaczenie ktore mamy sprawdzic
    public static double b = 0.06;
    public static boolean liczenie_bij_zLmax=true; //jesli na true progiem jest b/Lmax , jesli na false prog to b
    public static boolean wiazki_jednokierunkowe = true;
    public static double[][] AZ ={
            { 0.36429873, 0.5464481, 0.72859746, 0.0, 0.0, 0.36429873, },
            { 0.5464481, 0.81967217, 1.0928962, 0.0, 0.0, 0.5464481, },
            { 0.72859746, 1.0928962, 1.4571949, 0.0, 0.0, 0.72859746, },
            { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, },
            { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, },
            { 0.36429873, 0.5464481, 0.72859746, 0.0, 0.0, 0.36429873, },
    };


    //obl_Ilosc_Stanowisk
    public static float A=0.999f;
    public static float mi=1f/15f;
    public static float t=0.05f;
    public static float prawdopodobienstwo=0.1f;


    //Erlang i inne pierdoly
    public static float prog_erlang=0.01f;
    public static float A_erlang=1.5f;
    public static int N_erlang=5;



    public static float sum(float[] arr) {
        float tmp = 0;
        for (float v : arr)
            tmp += v;
        return tmp;
    }

    public static float[] mul(float[] src, float s) {
        float[] tmp = new float[src.length];
        for (int i = 0; i < src.length; i++) {
            tmp[i] = src[i] * s;
        }
        return tmp;
    }

    public static float[] mul(float[] src, float[] s) {
        float[] tmp = new float[src.length];
        for (int i = 0; i < src.length; i++) {
            tmp[i] = src[i] * s[i];
        }
        return tmp;
    }

    public static int[] getRoute(int from, int to) {

        Connection c = getConnection(connections_, from);
        for (Connection connection : c.connections) {
            if(connection.id==to)return new int[]{from,to};
        }

        ArrayList<Integer> arr = getRoute(c, new ArrayList<>(), to);
        int[] tmp = new int[arr.size()];
        for (int i = 0; i < arr.size(); i++) {
            tmp[i] = arr.get(i);
        }
        return tmp;
    }

    public static ArrayList<Integer> getRoute(Connection c, ArrayList<Integer> passedIDs, int to) {
        passedIDs.add(c.id);
        for (Connection connection : c.connections) {
            if (connection.id == to) {
                passedIDs.add(connection.id);
                return passedIDs;
            }
            if (!passedIDs.contains(connection.id)) {
                ArrayList<Integer> ps = new ArrayList<>(passedIDs);
                ArrayList<Integer> a = getRoute(connection, ps, to);
                if (a != null) return a;
            }
        }
        return null;
    }


    public static void main(String[] args) {
        if (typ == TypZadania.obl_AZ_Z_S_z_Roznymi_beta) {
            float[] Ai = mul(S, betaLista);
            float A = sum(Ai);
            Matrix AA = Matrix.diagonalFromArray(Ai);
            System.out.println("A");
            AA.print();
            System.out.println();
            Matrix wz = Matrix.WZ(sum(S) - 1, S);
            System.out.println("WZ");
            wz.print();
            System.out.println();
            System.out.println("AZ");
            Matrix matrix = Matrix.mult(AA, wz);
            matrix.print();
        } else if (typ == TypZadania.obl_AZ_Z_S) {
            float[] Ai = mul(S, beta);
            float A = sum(Ai);
            Matrix AA = Matrix.diagonalFromArray(Ai);
            System.out.println("A");
            AA.print();
            System.out.println();
            Matrix wz = Matrix.WZ(sum(S) - 1, S);
            System.out.println("WZ");
            wz.print();
            System.out.println();
            System.out.println("AZ");
            Matrix matrix = Matrix.mult(AA, wz);
            matrix.print();
        } else if (typ == TypZadania.obl_AZ_Z_Macierzy) {
            Matrix AA = Matrix.diagonalFromArray(Azm);
            Matrix filled = Matrix.fill(WZ_M);
            System.out.println("A");
            AA.print();
            System.out.println();
            System.out.println("WZ");
            filled.print();
            System.out.println();
            Matrix AZ = Matrix.mult(AA, filled);
            System.out.println("AZ");
            AZ.print();
        } else if (typ == TypZadania.obl_Zasoby_Wiazek) {
            connections();
            while(AWAIT){
                try {
                    Thread.sleep(100);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }

            Matrix matrix=new Matrix(AZ);
            System.out.println("AZ");
            matrix.print();
            polaczenie[0]--;
            polaczenie[1]--;
            int[] route = getRoute(polaczenie[0], polaczenie[1]);
            int[][] pairs = new int[route.length - 1][2];
            for (int i = 0; i < pairs.length; i++) {
                pairs[i][0] = route[i];
                pairs[i][1] = route[i + 1];
            }
            for (int[] pair : pairs) {
                System.out.println("nlz ( " + (pair[0] + 1) + ", " + (pair[1] + 1) + " )");
            }

            ArrayList<int[]> tmp = getAllCombinationsInList();

            float[] values = new float[pairs.length * 2];
            float[] Lmax = new float[2];
            int currentID = 0;
            for (int[] pair : pairs) {
                float tmpval = 0;
                System.out.println();
                System.out.print("azl( " + (pair[0] + 1) + ", " + (pair[1] + 1) + " ) = ");
                ArrayList<int[]> tmpp = new ArrayList<>();
                for (int[] ints : tmp) {
                    if (contains(pair, ints) && !arrayEquals(pair, ints)) {
                        tmpp.add(ints);
                        if (ints.length > Lmax[0]) Lmax[0] = ints.length;
                        tmpval += getValue(matrix, ints);
                        System.out.print(routeToString(ints) + " + ");
                    }
                }
                System.out.print(" = " + tmpval);
                values[currentID++] = tmpval;
                tmpp.clear();
                tmpval = 0;
                System.out.println();
                System.out.print("azl( " + (pair[1] + 1) + ", " + (pair[0] + 1) + " ) = ");
                for (int[] ints : tmp) {

                    if (containsSWAPPED(pair, ints) && !arrayEqualsSWAPPED(pair, ints)) {
                        tmpp.add(ints);
                        if (ints.length > Lmax[1]) Lmax[1] = ints.length;
                        tmpval += getValue(matrix, ints);
                        System.out.print(routeToString(ints) + " + ");
                    }
                }
                System.out.print(" = " + tmpval);
                values[currentID++] = tmpval;
            }
            boolean passed=false;
            for (float value : values) {
                if(value!=0)passed=true;
            }
            if(!passed){
                values[0]=getValue(matrix,new int[]{polaczenie[0],polaczenie[1]});
                values[1]=getValue(matrix,new int[]{polaczenie[1],polaczenie[0]});
            }

            Lmax[0]--;
            Lmax[1]--;
            System.out.println();
            System.out.println("Lmax[0]: "+Lmax[0]+" Lmax[1]: "+Lmax[1]);
                double bij1 = b / Lmax[0];
                double bij2 = b / Lmax[1];
                System.out.println();
                System.out.println(wiazki_jednokierunkowe?"Wiazki jednokieunkowe":"Wizaki dwukierunkowe");
                if(liczenie_bij_zLmax) {
                    System.out.println("bij1 = " + bij1);
                    System.out.println("bij2 = " + bij2);
                }
                int n = 1;
                for (int i = 0; i < pairs.length * 2; i++) {
                    if (i % 2 == 0) {
                        System.out.println();
                        System.out.println("Polaczenie: " + arrayToStringPlusOne(pairs[i / 2]));
                        float azlij = values[i]+(wiazki_jednokierunkowe?0:values[i+1]);
                        int N = 1;
                        double erlang = erlang(N, azlij);
                        if(pokaz_iteracje)System.out.println("Erlang N: " + N + " -> " + erlang);
                        while ((erlang = erlangRek(++N, azlij, erlang)) >(liczenie_bij_zLmax?(i >= pairs.length ? bij2 : bij1):b)) {
                            if (pokaz_iteracje) System.out.println("Erlang N: " + N + " -> " + erlang);
                        }
                        System.out.println(">>>Erlang( N: " + N + ", A: " + azlij + " ) = " + erlang);
                        System.out.println("N: " + N);
                    }
                }

        }else if(typ==TypZadania.obl_Ilosc_Stanowisk_Prawdopodobienstwo_Oczekiwania){
            int N=1;
            double erlang=erlang2Rek(1,A,0);
            if(pokaz_iteracje)System.out.println("P( tocz > "+t+") , N: " + N + " -> " + Pt(N,erlang));
            while ((Pt(N,erlang)) >prawdopodobienstwo) {
                erlang=erlang2Rek(++N,A,erlang);
                if (pokaz_iteracje)System.out.println("P( Tocz > "+t+") , N: " + N + " -> " + Pt(N,erlang));
            }
            System.out.println(">>>P( tocz > "+t+") , N: " + N + " -> " +  Pt(N,erlang));
            System.out.println("N: "+N);


        }else if(typ==TypZadania.obl_Pokaz_Wszystkie_Polaczenia_Przechodzace_Przez_Polaczenie){
            connections();
            while(AWAIT){
                try {
                    Thread.sleep(100);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }

            polaczenie[0]--;
            polaczenie[1]--;
            int[] route = getRoute(polaczenie[0], polaczenie[1]);
            int[][] pairs = new int[route.length - 1][2];
            for (int i = 0; i < pairs.length; i++) {
                pairs[i][0] = route[i];
                pairs[i][1] = route[i + 1];
            }
            for (int[] pair : pairs) {
                System.out.println("nlz ( " + (pair[0] + 1) + ", " + (pair[1] + 1) + " )");
            }

            ArrayList<int[]> tmp = getAllCombinationsInList();
            for (int[] pair : pairs) {
                System.out.println("Polaczenie: "+arrayToStringPlusOne(pair));
                for (int[] ints : tmp) {
                    if(contains(pair,route)){
                        System.out.println(arrayToStringPlusOne(ints));
                    }
                }
            }

        }else if(typ==TypZadania.obl_N_dla_progu_Erlang_2){
            int N = 1;
            double erlang = erlang2Rek(N, A_erlang);
            if(pokaz_iteracje)System.out.println("Erlang2 N: " + N + " -> " + erlang);
            while ((erlang = erlang2Rek(++N, A_erlang, erlang)) >prog_erlang) {
                if (pokaz_iteracje) System.out.println("Erlang2 N: " + N + " -> " + erlang);
            }
            System.out.println(">>>Erlang2( N: " + N + ", A: " + A_erlang + " ) = " + erlang);
            System.out.println("N: " + N);
        }else if(typ==TypZadania.obl_N_dla_progu_Erlang){
            int N = 1;
            double erlang = erlang(N, A_erlang);
            if(pokaz_iteracje)System.out.println("Erlang N: " + N + " -> " + erlang);
            while ((erlang = erlangRek(++N, A_erlang, erlang)) >prog_erlang) {
                if (pokaz_iteracje) System.out.println("Erlang N: " + N + " -> " + erlang);
            }
            System.out.println(">>>Erlang( N: " + N + ", A: " + A_erlang + " ) = " + erlang);
            System.out.println("N: " + N);
        }else if(typ==TypZadania.obl_Erlang_2_dla_N_A){
            System.out.println(">>>Erlang2( N: " + N_erlang + ", A: " + A_erlang + " ) = " + erlang2Rek(N_erlang,A_erlang));
        }else if(typ==TypZadania.obl_Erlang_dla_N_A){
            System.out.println(">>>Erlang( N: " + N_erlang + ", A: " + A_erlang + " ) = " + erlangRek(N_erlang,A_erlang));
        }

    }


    public static double Pt(int N,double erlang){
        return erlang*Math.pow(Math.E,-mi*(N-A)*t);
    }

    public static float getValue(Matrix mat, int[] route) {
        int i = route[0];
        int j = route[route.length - 1];
        return mat.val[i][j];
    }


    public static HashMap<Integer, ArrayList<int[]>> getAllCombinations() {
        HashMap<Integer, ArrayList<int[]>> tmp = new HashMap<Integer, ArrayList<int[]>>();
        for (int i = 0; i < connections_.size(); i++) {
            ArrayList<int[]> p = new ArrayList<>();
            for (int i1 = 0; i1 < connections_.size(); i1++) {
                if (i != i1) {
                    int[] route = getRoute(i, i1);
                    p.add(route);
                }
            }
            tmp.put(i, p);
        }
        return tmp;
    }

    public static ArrayList<int[]> getAllCombinationsInList() {
        ArrayList<int[]> tmp = new ArrayList<int[]>();
        for (int i = 0; i < connections_.size(); i++) {
            ArrayList<int[]> p = new ArrayList<>();
            for (int i1 = 0; i1 < connections_.size(); i1++) {
                if (i != i1) {
                    int[] route = getRoute(i, i1);
                    p.add(route);
                }
            }
            tmp.addAll(p);
        }
        return tmp;
    }

    static int suffixlen = 3, prefixlen = 5;

    public static void check(float f) {
        String[] tmp = ("" + f).replace(".", ";").split(";");
        String prefix = tmp[0];
        if (prefix.length() > prefixlen) {
            prefixlen = prefix.length();
        }
    }

    public static String routeToString(int[] route) {
        return "az( " + (route[0] + 1) + ", " + (route[route.length - 1] + 1) + " )";
    }

    public static boolean contains(int[] pair, int[] route) {
        for (int i = 0; i < route.length - 1; i++) {
            if (route[i] == pair[0] && route[i + 1] == pair[1]) return true;
        }
        return false;
    }

    public static boolean containsSWAPPED(int[] pair, int[] route) {
        for (int i = 0; i < route.length - 1; i++) {
            if (route[i] == pair[1] && route[i + 1] == pair[0]) return true;
        }
        return false;
    }


    public static String parse(float f) {
        String[] tmp = ("" + f).replace(".", ";").split(";");
        String prefix = "", suffix = tmp[1];
        if (tmp[0].length() > prefixlen) {
            System.out.println();
            System.out.println("======= BLAD COS SIE ROZJEBALO I CHUJ =========");
            System.out.println(f + ", " + tmp[0] + " " + prefixlen);
            System.out.println("Jakies dymy gosciu sie dzieja");
            System.exit(-1);
        }
        for (int i = 0; i < prefixlen - tmp[0].length(); i++) {
            prefix += " ";
        }
        prefix += tmp[0];
        if (suffix.length() >= suffixlen) {
            suffix = suffix.substring(0, suffixlen);
        } else {
            for (int i = 0; i <= suffixlen - suffix.length(); i++) {
                suffix += "0";
            }
        }
        return prefix + "." + suffix;
    }


    static class Matrix {
        float[][] val;
        int width, height;

        public Matrix(float[][] val) {
            this.val = val;
            this.height = val.length;
            this.width = val[0].length;
        }
        public Matrix(double[][] val) {
            this.val = doubleToFloatDoubleArray(val);
            this.height = val.length;
            this.width = val[0].length;
        }
        public static float multiplyMatricesCell(float[][] A, float[][] B, int row, int col) {
            float cell = 0;
            for (int i = 0; i < B.length; i++) {
                cell += A[row][i] * B[i][col];
            }
            return cell;
        }

        public static Matrix mult(Matrix A, Matrix B) {
            float[][] result = new float[A.height][B.width];
            for (int row = 0; row < result.length; row++) {
                for (int col = 0; col < result[row].length; col++) {
                    result[row][col] = multiplyMatricesCell(A.val, B.val, row, col);
                }
            }

            return new Matrix(result);


        }

        public static Matrix diagonal(int n) {
            float[][] m = new float[n][n];

            for (int i = 0; i < n; i++) m[i][i] = 1;

            return new Matrix(m);
        }

        public static Matrix AZ(float a, float[] ai) {
            float[][] m = new float[ai.length][ai.length];

            for (int i = 0; i < ai.length; i++) {
                for (int j = 0; j < ai.length; j++) {
                    m[i][j] = ai[i] * ai[j] / a;
                }
            }
            return new Matrix(m);
        }

        public static Matrix diagonalFromArray(float[] arr) {
            Matrix matrix = diagonal(arr.length);
            for (int i = 0; i < matrix.height; i++) {
                matrix.val[i] = mul(matrix.getRow(i), arr[i]);
            }
            return matrix;
        }

        public static Matrix WZ(float a, float[] s) {
            float[][] m = new float[s.length][s.length];

            for (int i = 0; i < s.length; i++) {
                for (int j = 0; j < s.length; j++) {
                    m[i][j] = s[j] / a;
                }
            }
            return new Matrix(m);
        }

        public static Matrix fill(double[][] aa) {
            for (int i = 0; i < aa.length; i++) {
                double[] row = aa[i];
                float tmp = 0;
                int id = -1;
                for (int i1 = 0; i1 < row.length; i1++) {
                    if (Math.abs(X - row[i1]) < 0.0001d) {
                        id = i1;
                    } else {
                        tmp += row[i1];
                    }

                }

                if (id != -1) {
                    if (tmp == 0) {
                        for (int i1 = 0; i1 < row.length; i1++) {
                            row[i1] = 0;
                        }
                    } else {
                        row[id] = 1f - tmp;
                    }
                }

                aa[i] = row;
            }
            return new Matrix(doubleToFloatDoubleArray(aa));
        }

        float[] getColumn(int id) {
            float[] h = new float[height];
            for (int j = 0; j < height; j++) {
                h[j] = val[j][id];
            }
            return h;
        }

        float[] getRow(int id) {
            return val[id];
        }


        public Matrix transpose() {

            float[][] v = new float[width][height];
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    v[x][y] = val[y][x];
                }
            }


            return new Matrix(v);
        }

        public static Matrix buildFromColumns(float[][] columns) {
            return new Matrix(columns).transpose();
        }


        public void print() {
            if(!LadneMacierze){
                printToCopy();
                return;
            }
            prefixlen = 0;
            System.out.println(height + " x " + width);
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    check(val[x][y]);
                }
            }
            for (int y = 0; y < height; y++) {
                System.out.print("| ");
                for (int x = 0; x < width; x++) {
                    System.out.print(parse(val[y][x]) + ", ");
                }
                System.out.println("|");
            }
        }

        public void printToCopy() {
            System.out.println("{");
            for (int y = 0; y < height; y++) {
                System.out.print("{ ");
                for (int x = 0; x < width; x++) {
                    System.out.print(val[y][x] + ", ");
                }
                System.out.println("},");
            }
            System.out.println("};");
        }
    }


    public static String arrayToStringPlusOne(int[] w) {
        String str = "[ ";
        for (int i = 0; i < w.length; i++) {
            str += ((w[i] + 1) + ", ");
        }
        str += "]";
        return str;
    }

    public static float[][] doubleToFloatDoubleArray(double[][] v) {
        float[][] a = new float[v.length][v[0].length];
        for (int i = 0; i < v.length; i++) {
            for (int i1 = 0; i1 < v[i].length; i1++) {
                a[i][i1] = ((float) v[i][i1]);
            }
        }
        return a;
    }

    public static Connection getConnection(ArrayList<Connection> connections, int id) {
        for (Connection connection : connections) {
            if (connection.id == id) return connection;
        }
        return null;
    }


    public static boolean arrayEquals(int[] a, int[] b) {
        if (a.length != b.length) return false;
        for (int i = 0; i < a.length; i++) {
            if (a[i] != b[i]) return false;
        }
        return true;
    }

    public static boolean arrayEqualsSWAPPED(int[] a, int[] b) {
        if (a.length != b.length) return false;
        for (int i = 0; i < a.length; i++) {
            if (a[a.length - i - 1] != b[i]) return false;
        }
        return true;
    }


    public static class Connection {
        int id;
        ArrayList<Connection> connections;

        public Connection(int id, ArrayList<Connection> connections) {
            this.id = id;
            this.connections = connections;
        }
    }

    public static double factorial(double num) {
        double tmp = 1;
        for (int i = 1; i <= num; i++) {
            tmp *= i;
        }
        return tmp;
    }

    public static double erlang(double N, double A) {
        double licznik = Math.pow(A, N) / factorial(N);
        double mianownik = 1;
        for (int i = 1; i <= N; i++) {
            mianownik += Math.pow(A, i) / factorial(i);
        }
        return licznik / mianownik;

    }

    public static double erlangRek(int N, double A, double erlangNM1) {
        if (N <= 1) return erlang(N, A);
        double licznik = A * erlangNM1;
        double mianownik = licznik + N;
        return licznik / mianownik;
    }

    public static double erlang2Rek(int N, double A, double erlangNM1) {
        if (N <= 1) return A;
        double licznik = A * (N-1-A)*erlangNM1;
        double mianownik = N*(N-1-A)+A*(1-erlangNM1);
        return licznik / mianownik;
    }

    public static double erlangRek(int N, double A) {
        if (N <= 1) return erlang(N, A);
        double erlangNM1=erlangRek(N-1,A);
        double licznik = A * erlangNM1;
        double mianownik = licznik + N;
        return licznik / mianownik;
    }

    public static double erlang2Rek(int N, double A) {
        if (N <= 1) return A;
        double erlangNM1=erlang2Rek(N-1,A);
        double licznik = A * (N-1-A)*erlangNM1;
        double mianownik = N*(N-1-A)+A*(1-erlangNM1);
        return licznik / mianownik;
    }


    //straszna sraka nizej
    public static void connections() {
        AWAIT = true;
        Thread t = new Thread() {
            ArrayList<Ab> points = new ArrayList<>();
            ArrayList<Con> connections = new ArrayList<>();

            class Ab {
                int id, x, y;

                public Ab(int id, int x, int y) {
                    this.id = id;
                    this.x = x;
                    this.y = y;
                }
            }

            class Con {
                Ab id1, id2;

                public Con(Ab id1, Ab id2) {
                    this.id1 = id1;
                    this.id2 = id2;
                }
            }

            public Ab get(int id) {
                for (Ab point : points) {
                    if (point.id == id) return point;
                }
                return null;
            }

            public float dst(Ab p, MouseEvent e) {
                float x = e.getX() - p.x - 25;
                float y = e.getY() - p.y - 25;
                return (float) Math.sqrt(x * x + y * y);
            }

            @Override
            public void run() {
                JFrame jFrame = new JFrame();
                jFrame.setSize(800, 600);
                int mx = 0;
                int my = 0;
                final int[] currentID = {0};
                final boolean[] conn = {false};
                final Ab[] current = {null};
                JPanel panel = new JPanel() {
                    @Override
                    public void paint(Graphics g) {
                        Graphics2D gg = (Graphics2D) g;
                        gg.clearRect(0, 0, 800, 600);
                        gg.drawString(!conn[0] ? "Wiazki cyz tam inne gowno" : "Polaczenia gosciu", 20, 20);
                        gg.drawString(!conn[0] ? "Wcisnij enter jak rozlozysz wezly" : "wcisnij enter zeby obliczyc wszystko", 20, 40);
                        gg.drawString("r - Resetuje wszystko, c - czysci " + (!conn[0] ? "abonentow" : "polaczenia"), 20, 60);
                        for (Ab point : points) {
                            gg.drawOval(point.x, point.y, 50, 50);
                            gg.drawString((point.id + 1) + "", point.x + 25, point.y + 25);
                        }
                        for (Con con : connections) {
                            gg.drawLine(con.id1.x + 25, con.id1.y + 25, con.id2.x + 25, con.id2.y + 25);
                        }
                    }
                };
                panel.addMouseListener(new MouseAdapter() {
                    @Override
                    public void mouseReleased(MouseEvent e) {
                        if (conn[0] && current[0] != null) {
                            Ab tmp = null;
                            for (Ab p : points) {
                                if (dst(p, e) < 25) {
                                    tmp = p;
                                }
                            }
                            if (tmp != null) {
                                if (tmp != current[0]) {
                                    connections.add(new Con(tmp, current[0]));
                                }
                                current[0] = null;
                                panel.repaint();
                            }
                        } else if (!conn[0]) {
                            points.add(new Ab(currentID[0]++, e.getX(), e.getY()));
                            panel.repaint();
                        }
                        super.mouseReleased(e);
                    }

                    @Override
                    public void mousePressed(MouseEvent e) {
                        if (conn[0] && current[0] == null) {
                            for (Ab p : points) {
                                if (dst(p, e) < 25) {
                                    current[0] = p;

                                }
                            }

                        }
                        super.mousePressed(e);
                    }
                });
                jFrame.addKeyListener(new KeyAdapter() {
                    @Override
                    public void keyPressed(KeyEvent e) {
                        if (e.getKeyCode() == KeyEvent.VK_ENTER) {
                            if (conn[0]) {
                                ArrayList<Connection> cn = new ArrayList<>();
                                for (Ab point : points) {
                                    cn.add(new Connection(point.id, new ArrayList<>()));
                                }
                                for (Con c : connections) {
                                    Connection c1 = getConnection(cn, c.id1.id);
                                    Connection c2 = getConnection(cn, c.id2.id);
                                    c1.connections.add(c2);
                                    c2.connections.add(c1);
                                }

                                connections_ = cn;
                                AWAIT = false;
                                jFrame.setVisible(false);
                                jFrame.dispose();
                            }
                            if (!conn[0]) conn[0] = true;
                            panel.repaint();
                        }
                        if (e.getKeyCode() == KeyEvent.VK_C) {
                            currentID[0] = 0;
                            if (conn[0]) {
                                connections.clear();
                            } else {
                                connections.clear();
                                points.clear();
                            }
                            panel.repaint();
                        }
                        if (e.getKeyCode() == KeyEvent.VK_R) {
                            connections.clear();
                            points.clear();
                            conn[0] = false;
                            currentID[0] = 0;
                            panel.repaint();
                        }
                        super.keyReleased(e);
                    }

                    @Override
                    public void keyReleased(KeyEvent e) {

                    }
                });


                jFrame.add(panel);
                jFrame.setVisible(true);
                jFrame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
            }
        };
        t.start();

    }


}