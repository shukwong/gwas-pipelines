version development

workflow move_outputs {
    input {
        File rclone_box_config
        File input_file
        String box_directory 
    }

    call move_file_to_box {
      input:
        File rclone_box_config
        File input_file
        String box_directory
    }

     output {
        String file_location_on_box = move_file_to_box.file_location_on_box
    }
}    

task move_file_to_box {
    
    File rclone_box_config
    File input

    String box_directory
    

    command <<<
		    rclone --config=~{rclone_box_config} copy ~{input} box:~{box_directory}
        rclone ls box:~{box_directory} > file_location_on_box.txt
    >>>    

	runtime {
		docker: "docker://rclone/rclone:1"
	}

    output {
	    String file_location_on_box = "file_location_on_box.txt"
    }
}

