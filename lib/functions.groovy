class functions {
    public static boolean checkFile(String path, String name, String extension) {

        def path_ = path.toString()
        def query = "find $path_/ -name $name*$extension"
        def answer = query.execute().text.count("\n") == 0 ? false : true

        return(answer)
    }
}